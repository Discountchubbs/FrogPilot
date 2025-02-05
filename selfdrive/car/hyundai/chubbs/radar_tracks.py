from collections.abc import Callable
import cereal.messaging as messaging
from openpilot.system import sentry
from openpilot.selfdrive.car.hyundai.values import HyundaiFlags, DBC
from openpilot.selfdrive.car.hyundai.radar_interface import RADAR_START_ADDR
from openpilot.selfdrive.car.hyundai.chubbs.fingerprint_utils import can_fingerprint
from openpilot.selfdrive.car.isotp_parallel_query import IsoTpParallelQuery


EXT_DIAG_REQUEST = b'\x10\x07'
EXT_DIAG_RESPONSE = b'\x50\x07'
WRITE_DATA_REQUEST = b'\x2e'
WRITE_DATA_RESPONSE = b'\x68'
RADAR_TRACKS_CONFIG = b"\x00\x00\x00\x01\x00\x01"

class RadarTrackController:
  def __init__(self, CP):
    self.CP = CP
    self.params = None

  def enable_radar_tracks(self, sendcan, logcan, next_can: Callable, bus=0, addr=0x7d0,
                         config_data_id=b'\x01\x42', timeout=0.1) -> bool:
    try:
      query = IsoTpParallelQuery(sendcan, next_can, bus, [addr], [EXT_DIAG_REQUEST], [EXT_DIAG_RESPONSE])
      for _, _ in query.get_data(timeout).items():
        query = IsoTpParallelQuery(sendcan, next_can, bus, [addr],
                                 [WRITE_DATA_REQUEST + config_data_id + RADAR_TRACKS_CONFIG],
                                 [WRITE_DATA_RESPONSE])
        query.get_data(0)
        return True
    except Exception:
      return False

  def log_fingerprint(self) -> None:
    # Log both dongle ID and car fingerprint to track which users enable radar tracks
    dongle_id = self.params.get("DongleId", encoding="utf-8")
    car_model = self.CP.carFingerprint
    fingerprint = [dongle_id, car_model]
    sentry.capture_message("Radar tracks enabled",
                          level="info",
                          extra={
                            "dongle_id": dongle_id,
                            "car_model": car_model,
                            "fingerprint": fingerprint
                          })

  def setup(self, params):
    self.params = params

  def initialize(self, params, logcan, sendcan) -> bool:
    if not self.params:
      self.setup(params)

    def get_next_can():
      return messaging.recv_one_retry(logcan)

    if self.CP.flags & HyundaiFlags.MANDO_RADAR and self.CP.radarUnavailable:
      if self.params.get_bool("HyundaiRadarTracksToggle"):
        self.CP.flags |= HyundaiFlags.ENABLE_RADAR_TRACKS
        if self.params.get_bool("HyundaiRadarTracks"):
          self.CP.radarUnavailable = False
          if self.enable_radar_tracks(sendcan, logcan, get_next_can):
            self.log_fingerprint()
            return True

    if self.CP.flags & HyundaiFlags.ENABLE_RADAR_TRACKS:
      _, fingerprint = can_fingerprint(get_next_can)
      radar_unavailable = RADAR_START_ADDR not in fingerprint[1] or self.CP.carFingerprint not in DBC

      radar_tracks = self.params.get_bool("HyundaiRadarTracks")
      radar_tracks_persistent = self.params.get_bool("HyundaiRadarTracksPersistent")

      self.params.put_bool_nonblocking("HyundaiRadarTracksConfirmed", radar_tracks)

      if not radar_tracks_persistent:
        self.params.put_bool_nonblocking("HyundaiRadarTracks", not radar_unavailable)
        self.params.put_bool_nonblocking("HyundaiRadarTracksPersistent", True)
        return not radar_unavailable

    return False
