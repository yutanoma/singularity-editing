#include "constrained_angles_cycle.h"

#include "soft_constrained_index_prescription.h"
#include "triangle_angles.h"

#include <array>
#include <igl/PI.h>

namespace rp {
namespace {
  double EPSILON = .0001;
  double OFFSET_LENGTH = 2.;
}

void add_brush_path_constraints(const Eigen::MatrixX3d &V, const Eigen::MatrixXi &EV,
                                const std::vector<std::vector<int>> &VE,
                                const std::vector<std::vector<double>> &VEAngles,
                                const std::vector<std::vector<EdgePath>> &brushPaths,
                                const Eigen::VectorXd &angleJumpPerEdge,
                                const Eigen::VectorXi &treeFathers,
                                const double &rotation,
                                std::vector<Eigen::Triplet<double>> & hardConstraintTriplets,
                                std::vector<double> &hardConstraintVals,
                                std::vector<Eigen::Triplet<double>> &softConstraintTriplets,
                                std::vector<double> &softConstraintVals) {
  for (int i = 0; i < brushPaths.size(); i++) {
    if (brushPaths[i].size() < 4) {
      continue;
    }

    int firstVid = -1;
    double firstEdgeVal = 0, firstEdgeAngle;

    // まずbrush内で同じ座標になるように
    // brushPathの端点は垂直とか定義できないので考慮しない。
    for (int j = 2; j < brushPaths[i].size() - 1; j++) {
      auto prev = brushPaths[i][j - 1];
      auto current = brushPaths[i][j];

      double prevVal = 0;
      double currentVal = 0;

      if (EV(prev.edgeId, 1) == EV(current.edgeId, 0)) {
        // 両方順方向
        prevVal = 1.0 - prev.ratio;
        currentVal = current.ratio;

        if (j == 2) {
          firstVid = EV(prev.edgeId, 0);
          firstEdgeVal = prev.ratio;
        }
      } else if (EV(prev.edgeId, 0) == EV(current.edgeId, 0)) {
        // prevが逆方向、currentが順方向
        prevVal = -prev.ratio;
        currentVal = current.ratio;

        if (j == 2) {
          firstVid = EV(prev.edgeId, 1);
          firstEdgeVal = prev.ratio - 1.0;
        }
      } else if (EV(prev.edgeId, 1) == EV(current.edgeId, 1)) {
        // pregが順方向、currentが逆方向
        prevVal = 1.0 - prev.ratio;
        currentVal = -1.0 + current.ratio;

        if (j == 2) {
          firstVid = EV(prev.edgeId, 0);
          firstEdgeVal = prev.ratio;
        }
      } else if (EV(prev.edgeId, 0) == EV(current.edgeId, 1)) {
        // prevもcurrentも逆方向
        prevVal = -prev.ratio;
        currentVal = -1.0 + current.ratio;

        if (j == 2) {
          firstVid = EV(prev.edgeId, 1);
          firstEdgeVal = prev.ratio - 1.0;
        }
      } else {
        assert("no path here");
      }

      // 【重要】進行方向左向きにベクトル場が向く！！！！
      // 辺の両端点は垂直に向いていると仮定して計算する
      std::array<double, 3> zeroForms = {.0, .0, .0};
      std::array<double, 2> oneForms = {.0, .0};
      std::array<int, 3> vids = {-1, -1, -1};
      std::array<double, 2> signs = {0, 0};

      for (int k = 0; k < 3; k++) {
        auto _prev = brushPaths[i][j - 2 + k];
        auto _next = brushPaths[i][j - 1 + k];

        int sharedVid = -1;
        std::array<int, 2> edgeIds = {};
        std::array<double, 2> lengths = {};
        for (int k = 0; k < 2; k++) {
          for (int l = 0; l < 2; l++) {
            if (EV(_prev.edgeId, k) == EV(_next.edgeId, l)) {
              sharedVid = EV(_prev.edgeId, k);
              double prevLength = (V.row(EV(_prev.edgeId, 1)) - V.row(EV(_prev.edgeId, 0))).norm();
              double nextLength = (V.row(EV(_next.edgeId, 1)) - V.row(EV(_next.edgeId, 0))).norm();

              if (k == 0) {
                lengths[0] = prevLength * _prev.ratio;
              } else if (k == 1) {
                lengths[0] = prevLength * (1 - _prev.ratio);
              } else {
                assert(false && "no ratio");
              }

              if (l == 0) {
                lengths[1] = nextLength * _next.ratio;
              } else if (l == 1) {
                lengths[1] = nextLength * (1 - _next.ratio);
              }

              edgeIds[0] = _prev.edgeId;
              edgeIds[1] = _next.edgeId;

              break;
            }
          }
        }
        assert(sharedVid != -1);

        Eigen::Vector3d prevp = (V.row(EV(_prev.edgeId, 0)) + (V.row(EV(_prev.edgeId, 1)) - V.row(EV(_prev.edgeId, 0))) * _prev.ratio).transpose();
        Eigen::Vector3d nextp = (V.row(EV(_next.edgeId, 0)) + (V.row(EV(_next.edgeId, 1)) - V.row(EV(_next.edgeId, 0))) * _next.ratio).transpose();
        double otherLen = (prevp - nextp).norm();

        double A, B, otherAngle;
        std::array<double, 2> angles = {.0, .0};
        rp::triangle_angles(lengths[0], lengths[1], otherLen, A, B, otherAngle);
        angles[0] = B;
        angles[1] = A;

        for (int l = 0; l < VE[sharedVid].size(); l++) {
          int length = VE[sharedVid].size();
          if (VE[sharedVid][l] == _prev.edgeId) {
            if (VE[sharedVid][(l + length - 1) % length] == _next.edgeId) {
              // next->prevの順
              double angle0 = VEAngles[sharedVid][(l + length - 1) % length];
              double angle1 = VEAngles[sharedVid][l];

              if (l == 0) {
                angle0 -= 2 * igl::PI;
              }

              assert((angle1 - angle0) >= 0);

              double angleRatio = (igl::PI / 2 - angles[1]) / otherAngle;
              zeroForms[k] = angle0 + (angle1 - angle0) * angleRatio;
            } else if (VE[sharedVid][(l + 1) % length] == _next.edgeId) {
              // prev->nextの順
              double angle0 = VEAngles[sharedVid][l];
              double angle1 = VEAngles[sharedVid][(l + 1) % length];

              if ((l + 1) % length == 0) {
                angle1 = 2 * igl::PI;
              }

              assert((angle1 - angle0) >= 0);

              double angleRatio = (igl::PI / 2 - angles[0]) / otherAngle;
              zeroForms[k] = igl::PI + (angle0 + (angle1 - angle0) * angleRatio);
            } else {
              assert(false && "next and prev not adjacent");
            }

            // std::cout << "[" << sharedVid << "] edgeAngles: {";
            // for (int m = 0; m < VE[sharedVid].size(); m++) {
            //   std::cout << VE[sharedVid][m] << " (" << VEAngles[sharedVid][m] * 180 / igl::PI << "), ";
            // }
            // std::cout << "} next: " << _next.edgeId << ", prev: " << _prev.edgeId << ", angle: " << zeroForms[k] * 180 / igl::PI << ", angles[0]: " << angles[0] * 180 / igl::PI << ", angles[1]: " << angles[1] * 180 / igl::PI << ", otherAngle: " << otherAngle * 180 / igl::PI << std::endl;

            vids[k] = sharedVid;

            break;
          }

          assert(l != VE[sharedVid].size() - 1);
        }

        if (k == 0) {
          if (EV(prev.edgeId, 0) == sharedVid) {
            signs[k] = 1;
          } else {
            signs[k] = -1;
          }
        } else if (k == 1) {
          if (EV(current.edgeId, 0) == sharedVid) {
            signs[k] = 1;
          } else {
            signs[k] = -1;
          }
        }
      }

      oneForms[0] = (zeroForms[1] - zeroForms[0]) * signs[0];
      oneForms[1] = (zeroForms[2] - zeroForms[1]) * signs[1];

      for (int k = 0; k < 2; k++) {
        if (vids[k] == vids[k+1]) {
          // 同じときは1-formはそのまま。
          oneForms[k] = k == 0 ? angleJumpPerEdge(prev.edgeId) : angleJumpPerEdge(current.edgeId);
        }
      }

      for (int k = 0; k < 2; k++) {
        double val = oneForms[k];
        double actual = k == 0 ? angleJumpPerEdge(prev.edgeId) : angleJumpPerEdge(current.edgeId);
        if (val - actual > igl::PI) {
          while (val - actual > igl::PI) {
            val -= 2 * igl::PI;
          }
        } else if (actual - val > igl::PI) {
          while (actual - val > igl::PI) {
            val += 2 * igl::PI;
          }
        }
        // std::cout << "original val: " << oneForms[k] << ", sign: " << signs[k] << ", val: " << val << ", actual: " << actual << ", diff: " << actual - val << std::endl;

        oneForms[k] = val;
      }

      if (prevVal != 0 || currentVal != 0) {
        hardConstraintTriplets.emplace_back(Eigen::Triplet<double>(hardConstraintVals.size(), prev.edgeId, prevVal));
        hardConstraintTriplets.emplace_back(Eigen::Triplet<double>(hardConstraintVals.size(), current.edgeId, currentVal));

        // std::cout << "prev: " << prevVal
        //           << ", next: " << currentVal
        //           << ", default angle jump0: " << angleJumpPerEdge(prev.edgeId) * prevVal
        //           << ", default angle jump1: " << angleJumpPerEdge(current.edgeId) * currentVal
        //           << ", new angle jump: " << oneForms[0] * prevVal + oneForms[1] * currentVal
        //           << ", diff: " << angleJumpPerEdge(prev.edgeId) * prevVal + angleJumpPerEdge(current.edgeId) * currentVal - (oneForms[0] * prevVal + oneForms[1] * currentVal) << std::endl;

        hardConstraintVals.emplace_back(oneForms[0] * prevVal + oneForms[1] * currentVal);
      }

      if (j == 2) {
        firstEdgeAngle = zeroForms[0] + oneForms[0] * firstEdgeVal;
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

    // std::cout << firstEdgeAngle << std::endl;

    softConstraintVals.emplace_back(firstEdgeAngle - rotation);
  }
}

void add_angle_constraints(const std::vector<int> &loop,
                           const Eigen::MatrixX3d &V,
                           const Eigen::MatrixXi &EV,
                           const std::vector<std::vector<int>> &VE,
                           const std::vector<std::vector<double>> &VEAngles,
                           const Eigen::VectorXi &treeFathers,
                           // signが正ならoutflow、負ならinflow
                           const int &sign,
                           double &globalRotation,
                           std::vector<Eigen::Triplet<double>> &hardConstraintMatrixTriplets,
                           std::vector<Eigen::Triplet<double>> &softConstraintMatrixTriplets,
                           std::vector<double> &angleJumpPerHardConstraint,
                           std::vector<double> &angleJumpPerSoftConstraint) {
  std::vector<int> edgeIds(loop.size()), diffs(loop.size());

  Eigen::VectorXd allAngles(loop.size());
  Eigen::VectorXi isBoundary(V.rows());
  isBoundary.setZero();

  for (int i = 0; i < loop.size(); i++) {
    int vid = loop[i];
    int nextVid = loop[(i + 1) % loop.size()];
    int edgeId = -1;

    isBoundary(vid) = 1;

    for (int j = 0; j < VE[vid].size(); j++) {
      int eid = VE[vid][j];
      if (EV(eid, 0) == vid) {
        edgeId = eid;
      } else if (EV(eid, 1) == vid) {
        edgeId = eid;
      }
    }

    assert(edgeId != -1);

    Eigen::Vector2d angles;
    double _sign = sign > 0 ? 1 : -1;

    for (int j = 0; j < 2; j++) {
      // VE[vid]は、0と最後が端になるようになっている
      int vid = EV(edgeId, j);
      double angle = VEAngles[vid][VEAngles[vid].size() - 1];
      if (angle < 0) {
        angle += 2 * igl::PI;
      }

      double vectorAngle = rp::round_pi(angle * _sign / 2);

      if (vid == 0) {
        globalRotation = vectorAngle;
      }

      // std::cout << "l323: [" << vid << "]: " << vectorAngle << ", " << _sign << ", " << angle << std::endl;

      if (vid == loop[i]) {
        allAngles(i) = vectorAngle;
      }
    
      angles(j) = vectorAngle;
    }

    double diff = - angles(0) + angles(1);
    edgeIds[i] = edgeId;
    diffs[i] = diff;
  }
  // std::cout << std::endl;

  Eigen::VectorXi isThinVertex(loop.size());
  isThinVertex.setZero();

  std::vector<int> treeLoopVids = {}, treeLoopIndices = {};
  for (int i = 0; i < loop.size(); i++) {
    int vid = loop[i];
    int nextVid = loop[(i + 1) % loop.size()];
    int prevVid = loop[(i + loop.size() - 1) % loop.size()];

    if ((V.row(nextVid) - V.row(prevVid)).norm() < EPSILON) {
      // 端から距離がlength以内だったら含めない
      double distance = 0;
      int count = 1;
      isThinVertex(i) = 1;
      int plusIndex = -1, plusVid = -1, minusIndex = -1, minusVid = -1;
      while (distance < OFFSET_LENGTH && count < loop.size() / 4) {
        plusIndex = (i + count) % loop.size();
        plusVid = loop[plusIndex];
        int plusPrevVid = loop[(i + count - 1) % loop.size()];

        minusIndex = (i + loop.size() - count) % loop.size();
        minusVid = loop[minusIndex];
        distance += (V.row(plusVid) - V.row(plusPrevVid)).norm();

        if (distance > OFFSET_LENGTH) {
          break;
        }

        isThinVertex(plusIndex) = 1;
        isThinVertex(minusIndex) = 1;
        count++;
      }
      if (plusIndex != -1 && plusVid != -1 && minusIndex != -1 && minusVid != -1) {
        minusIndex = (minusIndex + loop.size() - 1) % loop.size();
        minusVid = loop[minusIndex];
        treeLoopVids = {plusVid, minusVid};
        treeLoopIndices = {plusIndex, minusIndex};
      }
    }
  }

  int constraintsNum = 0;

  for (int i = 0; i < loop.size(); i++) {
    // std::cout << isThinVertex(i) << std::endl;
    if (isThinVertex(i) || isThinVertex((i + 1) % loop.size())) {
      continue;
    }

    if (constraintsNum >= loop.size()) {
      // 一周すると行列のランクが落ちるので抜ける
      continue;
    }

    int edgeId = edgeIds[i];
    double diff = diffs[i];
    int row = angleJumpPerHardConstraint.size();

    angleJumpPerHardConstraint.emplace_back(diff);

    hardConstraintMatrixTriplets.emplace_back(Eigen::Triplet<double>(row, edgeId, 1));

    constraintsNum++;
  }

  // std::cout << "l391: " << treeLoopVids.size() << std::endl;

  for (int i = 0; i < treeLoopVids.size(); i++) {
    int vid = treeLoopVids[i];
    // このvidから木を辿ってf0までの道のりを取る
    int currentVid = vid, prevVid;
    int count = angleJumpPerSoftConstraint.size();

    bool f0InPath = true;
    std::vector<Eigen::Triplet<double>> triplets = {};

    // std::cout << "l396: " << vid << std::endl;
    while (treeFathers(currentVid) != -1) {
      int edgeId = treeFathers(currentVid);
      double sign = EV(edgeId, 0) == currentVid ? -1 : 1;

      triplets.emplace_back(Eigen::Triplet<double>(count, edgeId, sign));

      currentVid = EV(edgeId, 0) == currentVid ? EV(edgeId, 1) : EV(edgeId, 0);
      std::cout << currentVid << ", ";

      if (!isBoundary(currentVid)) {
        f0InPath = false;
      }
    }
    std::cout << std::endl;

    if (f0InPath) {
      continue;
    }

    for (auto t : triplets) {
      softConstraintMatrixTriplets.emplace_back(t);
    }

    double diff = allAngles(treeLoopIndices[i]);
    angleJumpPerSoftConstraint.emplace_back(diff);
  }
}

void constrained_angles_cycle(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                              const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                              const std::vector<std::vector<int>> &VE,
                              const std::vector<std::vector<double>> &VEAngles,
                              const Eigen::VectorXi &treeFathers,
                              const Eigen::VectorXd &angleJumps,
                              const std::vector<std::vector<int>> &boundaryLoops,
                              const std::vector<std::vector<int>> &outflowPathVids,
                              const std::vector<std::vector<int>> &inflowPathVids,
                              const std::vector<std::vector<EdgePath>> &brushPaths,
                              const std::vector<Eigen::Triplet<double>> &basisLoopsTriplets,
                              const Eigen::SparseMatrix<double> &basisCyclesMatrix,
                              const int &N, const Eigen::SparseMatrix<double> &weightMatrix,
                              // ベクトル場とparameterizationの角度差
                              const double &rotation,
                              double &globalRotation,
                              Eigen::VectorXd &angleDefects,
                              Eigen::VectorXd &cycleIndices,
                              Eigen::SparseMatrix<double> &matrix,
                              Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &solver) {
  std::vector<Eigen::Triplet<double>> hardConstraintMatrixTriplets;
  std::vector<Eigen::Triplet<double>> softConstraintMatrixTriplets;

  std::vector<double> angleJumpPerHardConstraint = {};
  std::vector<double> angleJumpPerSoftConstraint = {};

  Eigen::VectorXi innerEdges(EV.rows());
  for (int i = 0; i < innerEdges.rows(); i++)
    innerEdges(i) = i;

  // outflowPathVidの場合は0はじまり、inflowPathVidの場合はoutflowPathVids.size()オフセットした値が入っている
  Eigen::VectorXi vid2PathIds(V.rows());
  vid2PathIds.setConstant(-1);
  for (int i = 0; i < outflowPathVids.size(); i++) {
    for (int j = 0; j < outflowPathVids[i].size(); j++) {
      vid2PathIds(outflowPathVids[i][j]) = i;
      // std::cout << i << ": " << outflowPathVids[i][j] << std::endl;
    }
  }

  for (int i = 0; i < inflowPathVids.size(); i++) {
    for (int j = 0; j < inflowPathVids[i].size(); j++) {
      vid2PathIds(inflowPathVids[i][j]) = i + outflowPathVids.size();
      // std::cout << i + outflowPathVids.size() << ": " << inflowPathVids[i][j] << std::endl;
    }
  }

  // 各boudnaryLoopに関して、外側/内側を判定してconstraintをかける
  for (int i = 0; i < boundaryLoops.size(); i++) {
    if (boundaryLoops[i].size() == 0) {
      continue;
    }

    int pathId = -1;
    for (int j = 0; j < boundaryLoops[i].size(); j++) {
      if (vid2PathIds(boundaryLoops[i][0]) != -1) {
        pathId = vid2PathIds(boundaryLoops[i][0]);
        break;
      }
    }

    if (pathId == -1) {
      // デフォルトはinflowにする
      std::cout << "sign!: -1: " << i << ", " << outflowPathVids.size() << ", " << pathId << std::endl;
      add_angle_constraints(boundaryLoops[i], V, EV, VE, VEAngles, treeFathers, -1, globalRotation, hardConstraintMatrixTriplets, softConstraintMatrixTriplets, angleJumpPerHardConstraint, angleJumpPerSoftConstraint);
    } else if (0 <= pathId && pathId < outflowPathVids.size()) {
      std::cout << "sign: 1: " << i << ", " << outflowPathVids.size() << ", " << pathId << std::endl;
      add_angle_constraints(boundaryLoops[i], V, EV, VE, VEAngles, treeFathers, 1, globalRotation, hardConstraintMatrixTriplets, softConstraintMatrixTriplets, angleJumpPerHardConstraint, angleJumpPerSoftConstraint);
    } else if (outflowPathVids.size() <= pathId) {
      std::cout << "sign: -1: " << i << ", " << outflowPathVids.size() << ", " << pathId << std::endl;
      add_angle_constraints(boundaryLoops[i], V, EV, VE, VEAngles, treeFathers, -1, globalRotation, hardConstraintMatrixTriplets, softConstraintMatrixTriplets, angleJumpPerHardConstraint, angleJumpPerSoftConstraint);
    }
  }

  std::cout << brushPaths.size() << std::endl;

  add_brush_path_constraints(V, EV, VE, VEAngles, brushPaths, angleJumps, treeFathers, rotation, hardConstraintMatrixTriplets, angleJumpPerHardConstraint, softConstraintMatrixTriplets, angleJumpPerSoftConstraint);

  for (int i = 0; i < angleJumpPerSoftConstraint.size(); i++) {
    angleJumpPerSoftConstraint[i] -= globalRotation;
  }

  // hardConstraintMatrixTriplets = {};
  // angleJumpPerHardConstraint = {};

  // softConstraintMatrixTriplets = {};
  // angleJumpPerSoftConstraint = {};

  rp::soft_constrained_index_prescription(V, F, EV, EF, innerEdges, angleJumps, N, weightMatrix, basisLoopsTriplets, basisCyclesMatrix, hardConstraintMatrixTriplets, angleJumpPerHardConstraint, softConstraintMatrixTriplets, angleJumpPerSoftConstraint, matrix, solver, angleDefects, cycleIndices);
};
}
