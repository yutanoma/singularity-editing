#pragma once

#include "./types.h"
#include "../store/geometry.h"
#include "../dispatcher/indices.h"

namespace rp {
namespace action {
EdgePathStateSetter use_brush_path_setter(ParameterizationStateSetter &parameterize, PathSetter &pathSetter);

EdgePathBrushStateSetter use_brush_stroke_setter(PathSetter &pathSetter);
}
}

