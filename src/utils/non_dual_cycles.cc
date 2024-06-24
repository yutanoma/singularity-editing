#include "./non_dual_cycles.h"

#include "./tree_cotree_decomposition.h"

#include <vector>
#include <Eigen/SparseCholesky>
#include <igl/boundary_loop.h>
#include <directional/tree.h>

#include <iostream>

namespace rp {
namespace non_dual_cycles {

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
              const std::vector<std::vector<int>> &VE,
              const Eigen::MatrixXi &FE,
              const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
              const std::vector<std::vector<int>> &boundaryLoops,
              Eigen::SparseMatrix<double> & basisCycles,
              std::vector<Eigen::Triplet<double>> &basisCycleTriplets)
{
  int eulerChar = V.rows() - EV.rows() + F.rows();

  int numBaseLoops = FE.rows();
  int numBoundaryLoops = boundaryLoops.size();
  int numGenerators = 2 - numBoundaryLoops - eulerChar;

  assert(numGenerators >= 0);

  std::cout << "numBaseLoops: " << numBaseLoops << ", numBoundaryLoops: " << numBoundaryLoops << ", numGenerators: " << numGenerators << ", " << "total: " << numBaseLoops + numBoundaryLoops + numGenerators << std::endl;

  basisCycles = Eigen::SparseMatrix<double>(numBaseLoops + numBoundaryLoops + numGenerators, EV.rows());
  basisCycleTriplets = {};

  Eigen::VectorXi isVisited(F.rows());
  isVisited.setZero();

  // 1. 各faceの回りのループ
  for (int i = 0; i < F.rows(); i++)
  {
    for (int j = 0; j < 3; j++) {
      int vid = F(i, j), nextVid = F(i, (j + 1) % 3);

      double sign = 0;
      int edgeId = -1;

      for (int k = 0; k < 3; k++) {
        if (EV(FE(i, k), 0) == vid && EV(FE(i, k), 1) == nextVid) {
          sign = +1.;
          edgeId = FE(i, k);
          break;
        } else if (EV(FE(i, k), 1) == vid && EV(FE(i, k), 0) == nextVid) {
          sign = -1.;
          edgeId = FE(i, k);
          break;
        }
      }

      assert(sign != 0 && edgeId != -1);

      basisCycleTriplets.emplace_back(Eigen::Triplet<double>(i, edgeId, sign));
    }
  }

  // 2. 各boundaryの回りのループ
  for (int i = 0; i < numBoundaryLoops; i++)
  {
    // std::cout << "[" << i << "](" << boundaryLoops[i].size() << "): ";
    for (int j = 0; j < boundaryLoops[i].size(); j++)
    {
      int vid = boundaryLoops[i][j];
      int nextVid = boundaryLoops[i][(j + 1) % boundaryLoops[i].size()];

      // std::cout << vid << ", ";

      int edgeId = -1;
      double sign = 0;
      for (int k = 0; k < VE[vid].size(); k++)
      {
        int eid = VE[vid][k];
        if (EV(eid, 0) == vid && EV(eid, 1) == nextVid)
        {
          edgeId = eid;
          // boundaryLoopはつねに時計回り
          sign = -1;
          break;
        } else if (EV(eid, 1) == vid && EV(eid, 0) == nextVid) {
          edgeId = eid;
          // boundaryLoopはつねに時計回りだが、このindexは反時計回り
          sign = 1;
          break;
        }
      }

      assert(edgeId != -1 && sign != 0);

      basisCycleTriplets.emplace_back(Eigen::Triplet<double>(numBaseLoops + i, edgeId, sign));
    }
    // std::cout << std::endl;
  }

  // 3. generator loops
  if (numGenerators > 0) {
    rp::tree_cotree_decomposition(EF, EV, numBaseLoops + numBoundaryLoops, basisCycleTriplets);
  }

  basisCycles.setFromTriplets(basisCycleTriplets.begin(), basisCycleTriplets.end());
}

void process(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
             const std::vector<std::vector<int>> &VE,
             const Eigen::MatrixXi &FE,
             const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EV,
             Eigen::SparseMatrix<double>& basisCycles) {
  std::vector<Eigen::Triplet<double>> basisCycleTriplets = {};
  std::vector<std::vector<int>> boundaryLoops;
  igl::boundary_loop(F, boundaryLoops);

  process(V, F, VE, FE, EF, EV, boundaryLoops, basisCycles, basisCycleTriplets);
}
}
}
