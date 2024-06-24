// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "dual_cycles.h"

#include "parallel_transport.h"
#include "vector_utils.h"

const double ROTATION_PER_EDGE_CRITERIA = 100;
namespace rp {
inline void get_edge_list_from_faces(const std::vector<int>& facePath,
                                     const Eigen::MatrixXi& FE,
                                     const Eigen::MatrixXi& EF,
                                     std::vector<int>& edgeList,
                                     std::vector<int>& edgeSignList) {
  edgeList.resize(facePath.size() - 1);
  edgeSignList.resize(facePath.size() - 1);

  for (int i = 0; i < facePath.size() - 1; i++) {
    int prevFid = facePath[i];
    int postFid = facePath[i + 1];

    int edgeId = -1;
    int edgeSign = -1;

    get_edge_from_faces(prevFid, postFid, FE, EF, edgeId, edgeSign);

    edgeList[i] = edgeId;
    edgeSignList[i] = edgeSign;
  }
};

void dual_cycles(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                 const Eigen::MatrixXi& EV, const Eigen::MatrixXi& EF,
                 Eigen::MatrixXi& FE, const std::vector<bool>& hasConstraint,
                 const std::vector<double>& constraint,
                 const Eigen::MatrixX3d& B1,
                 Eigen::SparseMatrix<double>& basisCycles,
                 Eigen::VectorXd& cycleCurvature, Eigen::VectorXi& vertex2cycle,
                 Eigen::VectorXi& innerEdges, int& innerVerticesSize) {
  using namespace Eigen;
  using namespace std;
  int numV = F.maxCoeff() + 1;
  int eulerChar = numV - EV.rows() + F.rows();
  vertex2cycle.conservativeResize(V.rows());

  std::vector<std::vector<int>> boundaryLoops;

  igl::boundary_loop(F, boundaryLoops);
  int numBoundaries = boundaryLoops.size();
  int numGenerators = 2 - numBoundaries - eulerChar;

  vector<Triplet<double>> basisCycleTriplets(EV.rows() * 2);
  basisCycles.resize(numV + numBoundaries + numGenerators, EV.rows());

  int numConstraints = 0;

  // all 1-ring cycles, including boundaries
  for (int i = 0; i < EV.rows(); i++) {
    basisCycleTriplets.push_back(Triplet<double>(EV(i, 0), i, -1.0));
    basisCycleTriplets.push_back(Triplet<double>(EV(i, 1), i, 1.0));
  }

  // Creating boundary cycles by building a matrix the sums up boundary loops
  // and zeros out boundary vertex cycles - it will be multiplied from the left
  // to basisCyclesMat
  VectorXi isBoundary(V.rows());
  isBoundary.setZero();
  for (int i = 0; i < boundaryLoops.size(); i++)
    for (int j = 0; j < boundaryLoops[i].size(); j++)
      isBoundary(boundaryLoops[i][j]) = 1;

  VectorXi pureInnerEdgeMask = VectorXi::Constant(EV.rows(), 1);
  for (int i = 0; i < EV.rows(); i++)
    if ((isBoundary(EV(i, 0))) || (isBoundary(EV(i, 1))))
      pureInnerEdgeMask(i) = 0;

  int currGeneratorCycle = 0;
  int currBoundaryCycle = 0;

  bool hasAnyConstraint = false;
  for (int i = 0; i < hasConstraint.size(); i++) {
    if (hasConstraint[i]) {
      hasAnyConstraint = true;
      break;
    }
  }

  vector<double> transportedConstraints = {};

  std::vector<std::vector<int>> wholeSignList;
  std::vector<std::vector<int>> wholeEdgeList;

  if ((numGenerators != 0 || hasAnyConstraint) /*||(numBoundaries!=0)*/) {
    MatrixXi reducedEV(EV);
    for (int i = 1; i < reducedEV.rows(); i++)
      if (isBoundary(reducedEV(i, 0)) || isBoundary(reducedEV(i, 1)))
        reducedEV(i, 0) = -1;

    VectorXi primalTreeEdges, primalTreeFathers;
    VectorXi dualTreeEdges, dualTreeFathers;
    directional::tree(reducedEV, primalTreeEdges, primalTreeFathers);
    // creating a set of dual edges that do not cross edges in the primal tree
    VectorXi fullIndices = VectorXi::LinSpaced(EV.rows(), 0, EV.rows() - 1);
    VectorXi reducedEFIndices, inFullIndices;
    MatrixXi reducedEF;
    igl::setdiff(fullIndices, primalTreeEdges, reducedEFIndices, inFullIndices);
    VectorXi Two = VectorXi::LinSpaced(2, 0, 1);

    igl::slice(EF, reducedEFIndices, Two, reducedEF);
    directional::tree(reducedEF, dualTreeEdges, dualTreeFathers);
    // dual tree edges is the spanning tree T
    // converting dualTreeEdges from reducedEF to EF
    for (int i = 0; i < dualTreeEdges.size(); i++)
      dualTreeEdges(i) = inFullIndices(dualTreeEdges(i));

    for (int i = 0; i < dualTreeFathers.size(); i++)
      if (dualTreeFathers(i) != -1 && dualTreeFathers(i) != -2)
        dualTreeFathers(i) = inFullIndices(dualTreeFathers(i));

    if (numGenerators != 0) {
      // building tree co-tree based homological cycles
      // finding dual edge which are not in the tree, and following their faces
      // to the end
      VectorXi isinTree = VectorXi::Zero(EF.rows());
      for (int i = 0; i < dualTreeEdges.size(); i++) {
        isinTree(dualTreeEdges(i)) = 1;
      }
      for (int i = 0; i < primalTreeEdges.size(); i++) {
        isinTree(primalTreeEdges(i)) = 1;
      }

      for (int i = 0; i < isinTree.size(); i++) {
        if (isinTree(i)) continue;

        // std::cout<<"New Cycle"<<std::endl;
        // otherwise, follow both end faces to the root and this is the dual
        // cycle
        if (EF(i, 0) == -1 || EF(i, 1) == -1) continue;
        std::vector<Triplet<double>> candidateTriplets;
        // candidateTriplets.push_back(Triplet<double>(0, i, 1.0));
        Vector2i currLeaves;
        currLeaves << EF(i, 0), EF(i, 1);
        VectorXi visitedOnce = VectorXi::Zero(
            EF.rows());  // used to remove the tail from the LCA to the root
        bool isBoundaryCycle = true;
        for (int i = 0; i < 2; i++) {  // on leaves
          int currTreeEdge = -1;       // indexing within dualTreeEdges
          int currFace = currLeaves(i);
          currTreeEdge = dualTreeFathers(currFace);
          if (currTreeEdge == -2) {
            break;
          }

          while (currTreeEdge != -1) {
            // std::cout<<"currTreeEdge: "<<currTreeEdge<<"\n"<<std::endl;
            // determining orientation of current edge vs. face
            double sign =
                ((EF(currTreeEdge, 0) == currFace) != (i == 0) ? 1.0 : -1.0);
            visitedOnce(currTreeEdge) = 1 - visitedOnce(currTreeEdge);
            candidateTriplets.push_back(Triplet<double>(0, currTreeEdge, sign));
            currFace = (EF(currTreeEdge, 0) == currFace ? EF(currTreeEdge, 1)
                                                        : EF(currTreeEdge, 0));
            currTreeEdge = dualTreeFathers(currFace);
          }
        }

        // only putting in dual edges that are below the LCA
        for (int i = 0; i < candidateTriplets.size(); i++)
          if ((visitedOnce(candidateTriplets[i].col())) &&
              (pureInnerEdgeMask(candidateTriplets[i].col())))
            isBoundaryCycle = false;

        int currRow =
            (isBoundaryCycle ? numV + currBoundaryCycle
                             : numV + numBoundaries + currGeneratorCycle);
        (isBoundaryCycle ? currBoundaryCycle++ : currGeneratorCycle++);

        basisCycleTriplets.push_back(Triplet<double>(currRow, i, 1.0));
        for (size_t i = 0; i < candidateTriplets.size(); i++)
          if (visitedOnce(candidateTriplets[i].col())) {
            Triplet<double> trueTriplet(currRow, candidateTriplets[i].col(),
                                        candidateTriplets[i].value());
            basisCycleTriplets.push_back(trueTriplet);
          }
      }
    }

    if (hasAnyConstraint && hasConstraint[0]) {
      // boundaryをconstraintにかける
      // 深さ優先探索

      std::vector<std::vector<int>> treeNextNodes;
      create_tree_from_tEf(dualTreeFathers, EF, treeNextNodes);

      std::stack<int> stack;
      stack.push(0);

      std::vector<bool> visited(F.size(), false);

      int constrainedPathsNum = 0,
          startIndex = numV + numBoundaries + numGenerators;

      // 深さ優先探索
      while (!stack.empty()) {
        int currentVid = stack.top();

        if (!visited[currentVid]) {
          if (hasConstraint[currentVid] && currentVid != 0) {
            // 前のconstrained valueを探す
            std::stack<int> copiedStack = stack;
            copiedStack.pop();

            std::vector<int> pathToAncestor = {}, edgeList = {},
                             edgeSignList = {};
            pathToAncestor.emplace_back(currentVid);

            int constrainedAncestor = -1;
            double confirmedDiff = 0;
            while (constrainedAncestor == -1) {
              int topOfStack = copiedStack.top();

              pathToAncestor.emplace_back(topOfStack);

              int prevFid = pathToAncestor[pathToAncestor.size() - 2];
              int postFid = pathToAncestor[pathToAncestor.size() - 1];

              int edgeId, edgeSign;

              get_edge_from_faces(prevFid, postFid, FE, EF, edgeId, edgeSign);

              edgeList.emplace_back(edgeId);
              edgeSignList.emplace_back(edgeSign);

              if (hasConstraint[topOfStack]) {
                double transportedConstraint =
                    parallel_transport(V, F, EF, FE, EV, B1, pathToAncestor,
                                       constraint[currentVid]);

                double diff = constraint[topOfStack] - transportedConstraint;

                // std::cout << diff << std::endl;

                if (diff > igl::PI) {
                  while (true) {
                    diff = diff - 2 * igl::PI;
                    if (diff <= igl::PI) {
                      break;
                    }
                  }
                } else if (diff < -igl::PI) {
                  while (true) {
                    diff = diff + 2 * igl::PI;
                    if (diff >= -igl::PI) {
                      break;
                    }
                  }
                }

                if (std::abs(diff / edgeList.size()) <
                    ROTATION_PER_EDGE_CRITERIA) {
                  constrainedAncestor = topOfStack;
                  confirmedDiff = diff;
                  break;
                }
              }

              copiedStack.pop();

              if (copiedStack.size() == 0) {
                constrainedAncestor = 0;
                edgeList.clear();
                edgeSignList.clear();
              }
            }

            if (edgeList.size() > 0) {
              // edge_sign_listをAに格納

              for (int i = 0; i < edgeSignList.size(); i++) {
                basisCycleTriplets.emplace_back(
                    Eigen::Triplet<double>(constrainedPathsNum + startIndex,
                                           edgeList[i], edgeSignList[i]));
              }

              // for debug
              wholeSignList.emplace_back(edgeSignList);
              wholeEdgeList.emplace_back(edgeList);

              // transportedConstraintをbに格納
              transportedConstraints.emplace_back(-confirmedDiff);

              constrainedPathsNum++;
            }

            // std::cout << diff * 180 / igl::PI << std::endl;
          }

          visited[currentVid] = true;
        }

        int nextVid = -1;

        for (int i = 0; i < treeNextNodes[currentVid].size(); i++) {
          int newVid = treeNextNodes[currentVid][i];
          ;

          if (!visited[newVid]) {
            nextVid = newVid;
            break;
          }
        }

        if (nextVid == -1) {
          stack.pop();
        } else {
          stack.emplace(nextVid);
        }
      }

      numConstraints = constrainedPathsNum;
    }
  }

  int matrixSize = numV + numBoundaries + numGenerators + numConstraints;

  std::cout << matrixSize << std::endl;

  SparseMatrix<double> sumBoundaryLoops(matrixSize, matrixSize);
  vector<Triplet<double>> sumBoundaryLoopsTriplets;
  vector<int> innerVerticesList, innerEdgesList;
  VectorXi remainRows, remainColumns;

  for (int i = 0; i < numV; i++) {
    sumBoundaryLoopsTriplets.push_back(
        Triplet<double>(i, i, 1.0 - isBoundary[i]));
    if (!isBoundary(i)) {
      innerVerticesList.push_back(i);
      vertex2cycle(i) = innerVerticesList.size() - 1;
    }
  }

  for (int i = 0; i < EV.rows(); i++)
    if (!((isBoundary(EV(i, 0))) && (isBoundary(EV(i, 1)))))
      innerEdgesList.push_back(i);

  // summing up boundary loops
  for (int i = 0; i < boundaryLoops.size(); i++)
    for (int j = 0; j < boundaryLoops[i].size(); j++) {
      sumBoundaryLoopsTriplets.push_back(
          Triplet<double>(numV + i, boundaryLoops[i][j], 1.0));
      vertex2cycle(boundaryLoops[i][j]) = innerVerticesList.size() + i;
    }

  // just passing generators through;
  for (int i = numV + numBoundaries; i < matrixSize; i++)
    sumBoundaryLoopsTriplets.push_back(Triplet<double>(i, i, 1.0));

  sumBoundaryLoops.setFromTriplets(sumBoundaryLoopsTriplets.begin(),
                                   sumBoundaryLoopsTriplets.end());

  basisCycles.resize(matrixSize, EV.rows());

  basisCycles.setFromTriplets(basisCycleTriplets.begin(),
                              basisCycleTriplets.end());

  basisCycles = sumBoundaryLoops * basisCycles;

  // removing rows and columns
  remainRows.resize(innerVerticesList.size() + numBoundaries + numGenerators +
                    numConstraints);

  remainColumns.resize(innerEdgesList.size());
  for (int i = 0; i < innerVerticesList.size(); i++)
    remainRows(i) = innerVerticesList[i];

  for (int i = 0; i < numBoundaries + numGenerators + numConstraints; i++)
    remainRows(innerVerticesList.size() + i) = numV + i;

  for (int i = 0; i < innerEdgesList.size(); i++)
    remainColumns(i) = innerEdgesList[i];

  // creating slicing matrices
  std::vector<Triplet<double>> rowSliceTriplets, colSliceTriplets;
  for (int i = 0; i < remainRows.size(); i++)
    rowSliceTriplets.push_back(Triplet<double>(i, remainRows(i), 1.0));
  for (int i = 0; i < remainColumns.size(); i++)
    colSliceTriplets.push_back(Triplet<double>(remainColumns(i), i, 1.0));

  SparseMatrix<double> rowSliceMat(remainRows.rows(), basisCycles.rows());
  rowSliceMat.setFromTriplets(rowSliceTriplets.begin(), rowSliceTriplets.end());

  SparseMatrix<double> colSliceMat(basisCycles.cols(), remainColumns.rows());
  colSliceMat.setFromTriplets(colSliceTriplets.begin(), colSliceTriplets.end());

  basisCycles = rowSliceMat * basisCycles * colSliceMat;

  cycleCurvature = VectorXd::Zero(basisCycles.rows());

  std::cout << basisCycles.rows() << " " << numGenerators << " "
            << numBoundaries << " " << numConstraints;

  innerEdges.conservativeResize(innerEdgesList.size());
  for (int i = 0; i < innerEdgesList.size(); i++)
    innerEdges(i) = innerEdgesList[i];

  // Correct computation of cycle curvature by adding angles
  // getting corner angle sum
  VectorXd allAngles(3 * F.rows());
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      RowVector3d edgeVec12 = V.row(F(i, (j + 1) % 3)) - V.row(F(i, j));
      RowVector3d edgeVec13 = V.row(F(i, (j + 2) % 3)) - V.row(F(i, j));
      allAngles(3 * i + j) =
          acos(edgeVec12.normalized().dot(edgeVec13.normalized()));
    }
  }

  // for each cycle, summing up all its internal angles negatively  + either
  // 2*pi*|cycle| for internal cycles or pi*|cycle| for boundary cycles.
  VectorXi isBigCycle = VectorXi::Ones(
      basisCycles.rows());  // TODO: retain it rather then reverse-engineer...

  for (int i = 0; i < V.rows(); i++)  // inner cycles
    if (!isBoundary(i)) isBigCycle(vertex2cycle(i)) = 0;

  // getting the 4 corners of each edge to allocated later to cycles according
  // to the sign of the edge.
  vector<set<int>> cornerSets(basisCycles.rows());
  vector<set<int>> vertexSets(basisCycles.rows());
  MatrixXi edgeCorners(innerEdges.size(), 4);
  for (int i = 0; i < innerEdges.rows(); i++) {
    int inFace1 = 0;
    while (F(EF(innerEdges(i), 0), inFace1) != EV(innerEdges(i), 0))
      inFace1 = (inFace1 + 1) % 3;
    int inFace2 = 0;
    while (F(EF(innerEdges(i), 1), inFace2) != EV(innerEdges(i), 1))
      inFace2 = (inFace2 + 1) % 3;

    edgeCorners(i, 0) = EF(innerEdges(i), 0) * 3 + inFace1;
    edgeCorners(i, 1) = EF(innerEdges(i), 1) * 3 + (inFace2 + 1) % 3;
    edgeCorners(i, 2) = EF(innerEdges(i), 0) * 3 + (inFace1 + 1) % 3;
    edgeCorners(i, 3) = EF(innerEdges(i), 1) * 3 + inFace2;
  }

  for (int k = 0; k < basisCycles.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(basisCycles, k); it; ++it) {
      cornerSets[it.row()].insert(
          edgeCorners(it.col(), it.value() < 0 ? 0 : 2));
      cornerSets[it.row()].insert(
          edgeCorners(it.col(), it.value() < 0 ? 1 : 3));
      vertexSets[it.row()].insert(
          EV(innerEdges(it.col()), it.value() < 0 ? 0 : 1));
    }

  for (int i = 0; i < cornerSets.size(); i++) {
    if (isBigCycle(i))
      cycleCurvature(i) = igl::PI * (double)(vertexSets[i].size());
    else
      cycleCurvature(i) = 2.0 * igl::PI;
    for (set<int>::iterator si = cornerSets[i].begin();
         si != cornerSets[i].end(); si++)
      cycleCurvature(i) -= allAngles(*si);
    for (int j = 0; j < transportedConstraints.size(); j++) {
      cycleCurvature(innerVerticesList.size() + numBoundaries + numGenerators +
                     j) = -transportedConstraints[j];
    }
  }

  innerVerticesSize = innerVerticesList.size();

  // for (int i = 0; i < 20; i++) {
  //   auto row = basisCycles.row(innerVerticesList.size() + numBoundaries +
  //                              numGenerators + i);
  //   auto elInVec = cycleCurvature(innerVerticesList.size() + numBoundaries +
  //                                 numGenerators + i);

  //   std::cout << "--------start:" << i << std::endl;

  //   std::cout << " !?" << std::endl;

  //   for (int j = 0; j < row.size(); j++) {
  //     if (basisCycles.coeff(
  //             innerVerticesList.size() + numBoundaries + numGenerators + i,
  //             j) != 0) {
  //       std::cout << j << " " << innerEdges(j) << " "
  //                 << basisCycles.coeff(innerVerticesList.size() +
  //                                          numBoundaries + numGenerators + i,
  //                                      j)
  //                 << std::endl;
  //     }
  //   }

  //   std::cout << "in vec"
  //             << ": " << elInVec;

  //   std::cout << "--actual result" << std::endl;

  //   for (int j = 0; j < wholeEdgeList[i].size(); j++) {
  //     std::cout << wholeEdgeList[i][j] << " " << wholeSignList[i][j]
  //               << std::endl;
  //   }

  //   std::cout << "in vec: " << transportedConstraints[i] << std::endl;
  // }
}
}  // namespace rp
