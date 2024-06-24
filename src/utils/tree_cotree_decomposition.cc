#include "tree_cotree_decomposition.h"
#include <iostream>

namespace rp {
void tree_cotree_decomposition(const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF, const int &rowOffset, std::vector<Eigen::Triplet<double>> &basisCycleTriplets) {
  using namespace Eigen;

  int currGeneratorCycle = 0;
  int currBoundaryCycle = 0;

  VectorXi primalTreeEdges, primalTreeFathers;
  VectorXi dualTreeEdges, dualTreeFathers;
  directional::tree(EV, primalTreeEdges, primalTreeFathers);
  VectorXi fullIndices = VectorXi::LinSpaced(EV.rows(), 0, EV.rows() - 1);
  VectorXi reducedEFIndices, inFullIndices;
  MatrixXi reducedEF;
  igl::setdiff(fullIndices, primalTreeEdges, reducedEFIndices, inFullIndices);
  VectorXi Two = VectorXi::LinSpaced(2, 0, 1);
  
  igl::slice(EF, reducedEFIndices, Two, reducedEF);
  directional::tree(reducedEF, dualTreeEdges, dualTreeFathers);
  //converting dualTreeEdges from reducedEF to EF
  for (int i = 0; i < dualTreeEdges.size(); i++)
    dualTreeEdges(i) = inFullIndices(dualTreeEdges(i));
  
  for (int i = 0; i < dualTreeFathers.size(); i++)
    if (dualTreeFathers(i) != -1 && dualTreeFathers(i) != -2)
      dualTreeFathers(i) = inFullIndices(dualTreeFathers(i));

  //building tree co-tree based homological cycles
  //finding dual edge which are not in the tree, and following their faces to the end
  VectorXi isinTree = VectorXi::Zero(EF.rows());
  for (int i = 0; i < dualTreeEdges.size(); i++) {
    isinTree(dualTreeEdges(i)) = 1;
  }
  for (int i = 0; i < primalTreeEdges.size(); i++) {
    isinTree(primalTreeEdges(i)) = 1;
  }

  for (int i = 0; i < isinTree.size(); i++) {
    if (isinTree(i))
      continue;
    
    //std::cout<<"New Cycle"<<std::endl;
    //otherwise, follow both end faces to the root and this is the dual cycle
    if (EF(i, 0) == -1 || EF(i, 1) == -1)
      continue;
    if (EV(i, 0) == -1 || EV(i, 1) == -1)
      continue;
    std::vector<Triplet<double> > candidateTriplets;
    //candidateTriplets.push_back(Triplet<double>(0, i, 1.0));
    Vector2i currLeaves; currLeaves << EF(i, 0), EF(i, 1);
    VectorXi visitedOnce = VectorXi::Zero(EF.rows());  //used to remove the tail from the LCA to the root

    for (int i = 0; i < 2; i++) { //on leaves
      int currTreeEdge = -1;  //indexing within dualTreeEdges
      int currFace = currLeaves(i);
      currTreeEdge = dualTreeFathers(currFace);
      if (currTreeEdge == -2) {
        break;
      }

      while (currTreeEdge != -1) {
        //std::cout<<"currTreeEdge: "<<currTreeEdge<<"\n"<<std::endl;
        //determining orientation of current edge vs. face
        double sign = ((EF(currTreeEdge, 0) == currFace) != (i == 0) ? 1.0 : -1.0);
        visitedOnce(currTreeEdge) = 1 - visitedOnce(currTreeEdge);
        candidateTriplets.push_back(Triplet<double>(0, currTreeEdge, sign));
        currFace = (EF(currTreeEdge, 0) == currFace ? EF(currTreeEdge, 1) : EF(currTreeEdge, 0));
        currTreeEdge = dualTreeFathers(currFace);
      }
    }

    int currRow = rowOffset + currGeneratorCycle;
    currGeneratorCycle++;

    basisCycleTriplets.push_back(Triplet<double>(currRow, i, 1.0));
    for (size_t i = 0; i < candidateTriplets.size(); i++)
      if (visitedOnce(candidateTriplets[i].col())) {
        Triplet<double> trueTriplet(currRow, candidateTriplets[i].col(), candidateTriplets[i].value());
        basisCycleTriplets.push_back(trueTriplet);
      }
  }
}
}
