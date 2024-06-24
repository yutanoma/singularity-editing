#pragma once

#include "./types.h"
#include "../dispatcher/button_operations.h"

namespace rp {
namespace action {
  RescalingStateSetter use_rescaling_setter(AngleSetter &angleSetter);

  RescalingStateSetter use_rescaling_resetter(AngleSetter &angleSetter);
}
}
