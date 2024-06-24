#pragma once

#include "../utils/opengl/glfw/Viewer.h"
#include <Eigen/Core>
#include <imgui/imgui.h>
#include "../utils/opengl/glfw/imgui/ImGuiMenu.h"
#include "../utils/opengl/glfw/imgui/ImGuiHelpers.h"
#include "../store/geometry.h"
#include "../action/types.h"
#include "./types.h"

namespace rp {
namespace viewer {
using namespace rp::action;
using WindowCallback = std::function<void()>;

WindowCallback use_custom_window(rp::opengl::glfw::Viewer &viewer,
                                 rp::opengl::glfw::imgui::ImGuiMenu &menu,
                                 rp::viewer::SelectMode &selectMode,
                                 bool &showVectorField,
                                 std::function<void()> &resetInteractionMode,
                                 VectorFieldStateSetter &renderVectorfield,
                                 std::function<void(void)> &clearVectorfield,
                                 ParameterizationStateSetter &parametrize,
                                 RescalingStateSetter &computeRescalings,
                                 RescalingStateSetter &clearRescalings);
}
}
