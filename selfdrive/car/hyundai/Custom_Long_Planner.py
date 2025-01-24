from openpilot.common.numpy_fast import clip, interp
from openpilot.common.conversions import Conversions as CV
from openpilot.common.realtime import DT_CTRL
from openpilot.selfdrive.car.hyundai.values import HyundaiFlags, Buttons, CANFD_CAR, LEGACY_SAFETY_MODE_CAR, CarControllerParams
from openpilot.selfdrive.car.hyundai import hyundaicanfd, hyundaican
from openpilot.selfdrive.controls.lib.longcontrol import LongCtrlState
from openpilot.selfdrive.controls.lib.drive_helpers import HYUNDAI_V_CRUISE_MIN

from openpilot.selfdrive.car.hyundai.hkg_additions import JerkLimiter
from openpilot.selfdrive.frogpilot.controls.frogpilot_planner import FrogPilotPlanner
from openpilot.selfdrive.frogpilot.controls.lib.frogpilot_acceleration import get_max_allowed_accel

class CustomStockLongitudinal:
  def __init__(self):
    self.frogpilot_planner = FrogPilotPlanner()
    self.jerk_limiter = JerkLimiter()

    # Button event tracking
    self.prev_cruise_buttons = []
    self.cruise_buttons = []

    # Button control variables
    self.timer = 0
    self.button_type = 0
    self.button_count = 0
    self.t_interval = 7
    self.cruise_button = None
    self.last_button_frame = 0
    self.frame = 0

    # State tracking
    self.road_curvature_detected = False
    self.v_set_dis = 0
    self.target_speed = 0
    self.v_tsc_state = 0
    self.m_tsc_state = 0
    self.traffic_mode_active = False
    self.v_cruise_min = 0

  def get_cruise_buttons_status(self, CS):
    if not CS.out.cruiseState.enabled or CS.cruise_buttons[-1] != Buttons.NONE:
      self.timer = 40
    elif self.timer:
      self.timer -= 1
    else:
      return 1
    return 0

  def get_button_type(self, CS, target_speed):
    current_speed = round(CS.out.cruiseState.speed *
                         (CV.MS_TO_MPH if not CS.params_list.is_metric else CV.MS_TO_KPH))

    if self.button_type == 0:
      self.button_count = 0
      if target_speed > current_speed:
        self.button_type = 1
      elif target_speed < current_speed:
        self.button_type = 2
      return None

    elif self.button_type == 1:
      self.button_count += 1
      if target_speed <= current_speed or self.button_count > 5:
        self.button_count = 0
        self.button_type = 3
      return Buttons.RES_ACCEL

    elif self.button_type == 2:
      self.button_count += 1
      if target_speed >= self.v_set_dis or self.v_set_dis <= self.v_cruise_min:
        self.button_count = 0
        self.button_type = 3
      elif self.button_count > 5:
        self.button_count = 0
        self.button_type = 3
      return Buttons.SET_DECEL

    elif self.button_type == 3:
      self.button_count += 1
      if self.button_count > self.t_interval:
        self.button_type = 0
      return None

    return None

  def update_custom_stock_long(self, car_fingerprint, frame, CS, packer, CP=None):
    self.frame = frame
    self.target_speed = CS.out.cruiseState.speed
    self.v_cruise_min = HYUNDAI_V_CRUISE_MIN[CS.params_list.is_metric] * (CV.KPH_TO_MPH if not CS.params_list.is_metric else 1)


    if not CS.out.cruiseState.enabled or not self.get_cruise_buttons_status(CS):
      return []

    if self.frogpilot_planner.sm.updated['frogpilotPlan']:
      fp = self.frogpilot_planner.sm['frogpilotPlan']
      self.road_curvature_detected = fp.vtscControllingCurve
      self.traffic_mode_active = fp.trafficModeActive
      self.v_tsc_state = fp.vtscSpeed
      self.m_tsc_state = fp.mtscSpeed
      self.speed_limit = fp.slcSpeedLimit

    MS_CONVERT = CV.MS_TO_KPH if CS.params_list.is_metric else CV.MS_TO_MPH
    v_cruise_target = round(self.frogpilot_planner.v_cruise * MS_CONVERT)
    self.v_set_dis = round(CS.out.cruiseState.speed * MS_CONVERT)

    self.cruise_button = self.get_button_type(CS, v_cruise_target)

    if self.cruise_button is None:
        return []

    if (self.jerk_limiter.hkg_tuning and CP) or self.frogpilot_planner.sport_plus:
      if self.frogpilot_planner.sport_plus:
        accel = min(self.frogpilot_planner.max_desired_acceleration, get_max_allowed_accel(CS.out.vEgo))
      else:
        accel = min(self.frogpilot_planner.max_desired_acceleration, CarControllerParams.ACCEL_MAX)

      if accel < 0 or self.jerk_limiter.using_e2e:
        accel = self.jerk_limiter.calculate_limited_accel(accel, CS, LongCtrlState, interp, clip)

    return self.get_hyundai_messages(car_fingerprint, CS, packer, CP, v_cruise_target)

  def get_hyundai_messages(self, car_fingerprint, CS, packer, CP, v_cruise_target):
    """Create Hyundai-specific messages for longitudinal control"""
    if car_fingerprint in CANFD_CAR:
      if CP.flags & HyundaiFlags.CANFD_ALT_BUTTONS:
        return []

      if self.frame % 2 == 0:
        return [hyundaicanfd.create_buttons(packer, CP, CP.CAN,
                ((self.frame // 2) + 1) % 0x10, self.cruise_button)]

    else:
      if car_fingerprint in LEGACY_SAFETY_MODE_CAR:
        traffic_factor = 2 if self.traffic_mode_active else 1
        base_freq = 5 if abs(v_cruise_target - self.v_set_dis) <= 2 else 1
        send_freq = base_freq * traffic_factor

        if (self.frame - self.last_button_frame) * DT_CTRL > 0.1 * send_freq:
          messages = [hyundaican.create_clu11(packer, self.frame, CS.clu11, self.cruise_button)] * 25

          if (self.frame - self.last_button_frame) * DT_CTRL >= 0.15 * send_freq:
            self.last_button_frame = self.frame

          return messages
      else:
        if self.frame % 2 == 0:
          return [hyundaican.create_clu11(packer, (self.frame // 2) + 1, CS.clu11, self.cruise_button)] * 25

    return []
