#include "./window.h"

#include "./mouse_down.h"

namespace rp {
namespace viewer {
namespace {
  rp::viewer::SelectMode innerSelectMode = rp::viewer::SelectMode::undefined;
  bool innerShowVectorField = false;
}

void basic_setup(rp::opengl::glfw::imgui::ImGuiMenu &menu) {
  // Define next window position + size
  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 0), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
  ImGui::Begin(
      "Singularity Editing  ", nullptr,
      ImGuiWindowFlags_NoSavedSettings
  );

  // Expose the same variable directly ...
  ImGui::PushItemWidth(-80);
  ImGui::PopItemWidth();
}

void radiobutton_operations(rp::viewer::SelectMode &selectMode, std::function<void()> &resetInteractionMode) {
  int _selectMode = selectMode;

  ImGui::RadioButton("move singularities", &_selectMode, rp::viewer::SelectMode::exchangeIndices);
  ImGui::RadioButton("add singularity pairs", &_selectMode, rp::viewer::SelectMode::addSingularities);
  ImGui::RadioButton("remove singularities", &_selectMode, rp::viewer::SelectMode::annihilateIndices);
  // ImGui::RadioButton("add brush path", &_selectMode, rp::viewer::SelectMode::brushPath);
  // ImGui::RadioButton("add index prescription path", &_selectMode, rp::viewer::SelectMode::indexPrescriptionPath);

  switch (_selectMode) {
    case rp::viewer::SelectMode::exchangeIndices:
      selectMode = rp::viewer::SelectMode::exchangeIndices;
      break;
    case rp::viewer::SelectMode::annihilateIndices:
      selectMode = rp::viewer::SelectMode::annihilateIndices;
      break;
    case rp::viewer::SelectMode::addSingularities:
      selectMode = rp::viewer::SelectMode::addSingularities;
      break;
    // todo: make this work
    // case rp::viewer::SelectMode::brushPath:
    //   selectMode = rp::viewer::SelectMode::brushPath;
    //   break;
    // case rp::viewer::SelectMode::indexPrescriptionPath:
    //   selectMode = rp::viewer::SelectMode::indexPrescriptionPath;
    //   break;
    default:
      selectMode = rp::viewer::SelectMode::undefined;
      break;
  }

  if (innerSelectMode != selectMode) {
    resetInteractionMode();
    innerSelectMode = selectMode;
  }
}

void checkbox_operations(bool &showVectorField, VectorFieldStateSetter &renderVectorfield, std::function<void(void)> &clearVectorfield) {
  ImGui::Checkbox("show vector field", &showVectorField);

  if (showVectorField != innerShowVectorField) {
    if (showVectorField) {
      renderVectorfield();
    } else {
      clearVectorfield();
    }
    innerShowVectorField = showVectorField;
  }
}

void button_operations(RescalingStateSetter &computeRescalings,
                       RescalingStateSetter &clearRescalings,
                       ParameterizationStateSetter &parametrize) {
  if (ImGui::Button("compute rescaling", ImVec2()))
  {
    computeRescalings();
  }

  if (ImGui::Button("clear rescaling", ImVec2()))
  {
    clearRescalings();
  }

  if (ImGui::Button("parameterize", ImVec2())) {
    parametrize();
  }
}

WindowCallback use_custom_window(rp::opengl::glfw::Viewer &viewer, 
                                 rp::opengl::glfw::imgui::ImGuiMenu &menu,
                                 rp::viewer::SelectMode &selectMode,
                                 bool &showVectorField,
                                 std::function<void()> &resetInteractionMode,
                                 VectorFieldStateSetter &renderVectorfield,
                                 std::function<void(void)> &clearVectorfield,
                                 ParameterizationStateSetter &parametrize,
                                 RescalingStateSetter &computeRescalings,
                                 RescalingStateSetter &clearRescalings) {
  auto windowCallback = [&]()
  {
    basic_setup(menu);

    radiobutton_operations(selectMode, resetInteractionMode);

    checkbox_operations(showVectorField, renderVectorfield, clearVectorfield);

    button_operations(computeRescalings, clearRescalings, parametrize);

    ImGui::End();
  };

  return windowCallback;
};
};
};
