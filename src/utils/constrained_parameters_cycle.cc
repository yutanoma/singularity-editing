#include "constrained_parameters_cycle.h"

#include "soft_constrained_index_prescription.h"

#include <igl/PI.h>

namespace rp {
  void add_index_prescription_constraints(const Eigen::MatrixXi &EV,
                                          const std::vector<std::vector<EdgePath>> &indexPrescriptionPaths,
                                          const std::vector<double> indexPrescriptionPathIndices,
                                          std::vector<Eigen::Triplet<double>> &hardConstraintTriplets,
                                          std::vector<double> &hardConstraintVals) {
    for (int i = 0; i < indexPrescriptionPaths.size(); i++) {
      for (int j = 0; j < indexPrescriptionPaths[i].size(); j++) {
        // 前→後が正方向
        auto current = indexPrescriptionPaths[i][j];
        double val = 1.0;

        if (j > 0) {
          auto prev = indexPrescriptionPaths[i][j - 1];
          if (EV(current.edgeId, 0) == EV(prev.edgeId, 0) || EV(current.edgeId, 0) == EV(prev.edgeId, 1)) {
            val = 1.0;
          } else {
            val = -1.0;
          }
        } else if (j == 0) {
          auto next = indexPrescriptionPaths[i][j + 1];
          if (EV(current.edgeId, 0) == EV(next.edgeId, 0) || EV(current.edgeId, 0) == EV(next.edgeId, 1)) {
            val = -1.0;
          } else {
            val = 1.0;
          }
        }

        if (j == 0) {
          if (val > 0) {
            // 順方向
            val = 1.0 - current.ratio;
          } else {
            val = -current.ratio;
          }
        } else if (j == indexPrescriptionPaths[i].size() - 1) {
          if (val > 0) {
            val = current.ratio;
          } else {
            val = current.ratio - 1.0;
          }
        }

        hardConstraintTriplets.emplace_back(Eigen::Triplet<double>(hardConstraintVals.size(), current.edgeId, val));
      }

      hardConstraintVals.emplace_back(0);
    }
  }

  void add_brush_paths_constraints(const Eigen::MatrixXi &EV,
                                   const std::vector<std::vector<EdgePath>> &brushPaths,
                                   const double &levelsetVal,
                                   const Eigen::VectorXd &projectionPerEdge,
                                   const Eigen::VectorXi &treeFathers,
                                   std::vector<Eigen::Triplet<double>> & hardConstraintTriplets,
                                   std::vector<double> &hardConstraintVals,
                                   std::vector<Eigen::Triplet<double>> &softConstraintTriplets,
                                   std::vector<double> &softConstraintVals) {
    // std::cout << "brushpathsize: " << brushPaths.size() << std::endl;
    for (int i = 0; i < brushPaths.size(); i++) {
      // std::cout << "aaa " << brushPaths[i].size() << std::endl;
      if (brushPaths[i].size() <= 1) {
        continue;
      }

      int firstVid = -1;
      double firstEdgeVal = 0;

      // まずbrush内で同じ座標になるように
      for (int j = 1; j < brushPaths[i].size(); j++) {
        auto prev = brushPaths[i][j - 1];
        auto current = brushPaths[i][j];

        double prevVal = 0;
        double currentVal = 0;

        if (EV(prev.edgeId, 1) == EV(current.edgeId, 0)) {
          // 両方順方向
          prevVal = 1.0 - prev.ratio;
          currentVal = current.ratio;

          if (j == 1) {
            firstVid = EV(prev.edgeId, 0);
            firstEdgeVal = prev.ratio;
          }
        } else if (EV(prev.edgeId, 0) == EV(current.edgeId, 0)) {
          // prevが逆方向、currentが順方向
          prevVal = -prev.ratio;
          currentVal = current.ratio;

          if (j == 1) {
            firstVid = EV(prev.edgeId, 1);
            firstEdgeVal = prev.ratio - 1.0;
          }
        } else if (EV(prev.edgeId, 1) == EV(current.edgeId, 1)) {
          // pregが順方向、currentが逆方向
          prevVal = 1.0 - prev.ratio;
          currentVal = -1.0 + current.ratio;

          if (j == 1) {
            firstVid = EV(prev.edgeId, 0);
            firstEdgeVal = prev.ratio;
          }
        } else if (EV(prev.edgeId, 0) == EV(current.edgeId, 1)) {
          // prevもcurrentも逆方向
          prevVal = -prev.ratio;
          currentVal = -1.0 + current.ratio;

          if (j == 1) {
            firstVid = EV(prev.edgeId, 1);
            firstEdgeVal = prev.ratio - 1.0;
          }
        } else {
          assert("no path here");
        }

        if (prevVal != 0 || currentVal != 0) {
          hardConstraintTriplets.emplace_back(Eigen::Triplet<double>(hardConstraintVals.size(), prev.edgeId, prevVal));
          hardConstraintTriplets.emplace_back(Eigen::Triplet<double>(hardConstraintVals.size(), current.edgeId, currentVal));

          hardConstraintVals.emplace_back(0);
        }
      }

      assert(firstVid != -1);

      // 最初のedgeについて、最初のedge上の点までのtree上のパスを作る

      // このvidから木を辿ってf0までの道のりを取る
      int currentVid = firstVid, prevVid;
      int count = softConstraintVals.size();
      while (treeFathers(currentVid) != -1) {
        int edgeId = treeFathers(currentVid);
        double sign = EV(edgeId, 0) == currentVid ? -1 : 1;

        softConstraintTriplets.emplace_back(Eigen::Triplet<double>(count, edgeId, sign));

        currentVid = EV(edgeId, 0) == currentVid ? EV(edgeId, 1) : EV(edgeId, 0);
      }

      softConstraintTriplets.emplace_back(Eigen::Triplet<double>(count, brushPaths[i][0].edgeId, firstEdgeVal));

      softConstraintVals.emplace_back(levelsetVal);
    }
  }

  // パスの種類は
  // 1. indexPrescriptionPaths: ブラシで書いた場所のindexが指定した値になるようにするパス。パスの内部はhardConstraint
  // 2. brushPaths: ブラシで描いた場所がレベルセットの界面になるようにするパス。パス内部はhardConsrtaint、パスの端点とf0の間はsoftConstraint
  void constrained_parameters_cycle(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                                    const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                                    const Eigen::VectorXi &treeFathers,
                                    const Eigen::VectorXd &projectionPerEdge,
                                    const std::vector<Eigen::Triplet<double>> &basisLoopsTriplets,
                                    const Eigen::SparseMatrix<double> &basisCyclesMatrix,
                                    const std::vector<std::vector<EdgePath>> &indexPrescriptionPaths,
                                    const std::vector<double> indexPrescriptionPathIndices,
                                    const std::vector<std::vector<EdgePath>> &brushPaths,
                                    const int &N, const Eigen::SparseMatrix<double> &weightMatrix,
                                    const double &levelsetVal,
                                    double &globalRotation,
                                    Eigen::VectorXd &cycleIndices,
                                    Eigen::SparseMatrix<double> &matrix,
                                    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver) {
    std::vector<Eigen::Triplet<double>> hardConstraintMatrixTriplets;
    std::vector<Eigen::Triplet<double>> softConstraintMatrixTriplets;

    std::vector<double> valueJumpPerHardConstraint = {};
    std::vector<double> valueJumpPerSoftConstraint = {};

    Eigen::VectorXi innerEdges(EV.rows());
    for (int i = 0; i < innerEdges.rows(); i++)
      innerEdges(i) = i;

    add_index_prescription_constraints(EV, indexPrescriptionPaths, indexPrescriptionPathIndices, hardConstraintMatrixTriplets, valueJumpPerHardConstraint);
    add_brush_paths_constraints(EV, brushPaths, levelsetVal, projectionPerEdge, treeFathers, hardConstraintMatrixTriplets, valueJumpPerHardConstraint, softConstraintMatrixTriplets, valueJumpPerSoftConstraint);

    for (int i = 0; i < valueJumpPerSoftConstraint.size(); i++) {
      valueJumpPerSoftConstraint[i] -= globalRotation;
    }

    Eigen::VectorXd lengthDefects = basisCyclesMatrix * projectionPerEdge;

    rp::soft_constrained_index_prescription(V, F, EV, EF, innerEdges, projectionPerEdge, N, weightMatrix, basisLoopsTriplets, basisCyclesMatrix, hardConstraintMatrixTriplets, valueJumpPerHardConstraint, softConstraintMatrixTriplets, valueJumpPerSoftConstraint, matrix, solver, lengthDefects, cycleIndices);
  }
}
