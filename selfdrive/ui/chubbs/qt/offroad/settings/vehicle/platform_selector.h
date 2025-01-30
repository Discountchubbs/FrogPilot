/**
 * Copyright (c) 2021-, Haibin Wen, chubbs, and a number of other contributors.
 *
 * This file is part of chubbs and is licensed under the MIT License.
 * See the LICENSE.md file in the root directory for more details.
 */

#pragma once

#include "selfdrive/ui/chubbs/qt/offroad/settings/settings.h"

class PlatformSelector : public ButtonControl {
  Q_OBJECT

public:
  PlatformSelector();
  QVariant getPlatformBundle(const QString &key);

public slots:
  void refresh(bool _offroad);

private:
  void searchPlatforms(const QString &query);
  QMap<QString, QVariantMap> loadPlatformList();

  Params params;
  bool offroad;
};
