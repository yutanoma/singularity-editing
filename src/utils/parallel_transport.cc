#include "parallel_transport.h"

#include "vector_utils.h"
#include <iostream>

namespace rp {
double parallel_transport(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                          const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                          const Eigen::MatrixXi &EV, const Eigen::MatrixX3d &B1,
                          const std::vector<int> &face_path,
                          const double &initial_angle) {
  std::vector<int> edgeIds = {}, edgeSigns = {};

  return parallel_transport(V, F, EF, FE, EV, B1, face_path, initial_angle, edgeIds, edgeSigns);
};

double edges_parallel_transport(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                                const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                                const Eigen::MatrixXi &EV, const Eigen::MatrixX3d &B1,
                                const std::vector<int> &edge_path, const std::vector<int> &edge_signs,
                                const double &initial_angle) {
  double transported_angle = initial_angle;

  for (int i = 0; i < edge_path.size(); i++) {
    int edge_id = edge_path[i];

    int prev_fid = edge_signs[i] > 0 ? EF(edge_id, 0) : EF(edge_id, 1);
    int post_fid = edge_signs[i] > 0 ? EF(edge_id, 1) : EF(edge_id, 0);

    Eigen::Vector3d edgeVec =
        (V.row(EV(edge_id, 1)).transpose() - V.row(EV(edge_id, 0)).transpose())
            .normalized();
    Eigen::Vector3d invertedEdgeVec =
        (V.row(EV(edge_id, 0)).transpose() - V.row(EV(edge_id, 1)).transpose())
            .normalized();

    Eigen::Vector3d prevNormal = face_normal(V, F, prev_fid);
    Eigen::Vector3d postNormal = face_normal(V, F, post_fid);

    double prev_angle =
        get_angle(B1.row(prev_fid).transpose(), edgeVec, prevNormal);
    double post_angle =
        get_angle(B1.row(post_fid).transpose(), edgeVec, postNormal);

    double edge_jump = -prev_angle + post_angle;

    double inverted_prev_angle =
        get_angle(B1.row(prev_fid).transpose(), invertedEdgeVec, prevNormal);
    double inverted_post_angle =
        get_angle(B1.row(post_fid).transpose(), invertedEdgeVec, postNormal);

    double inverted_edge_jump = -inverted_prev_angle + inverted_post_angle;

    edge_jump = std::abs(edge_jump) > std::abs(inverted_edge_jump) ? inverted_edge_jump : edge_jump;

    transported_angle = transported_angle - prev_angle + post_angle;
  }

  return transported_angle;
}

double parallel_transport(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
                          const Eigen::MatrixXi &EF, const Eigen::MatrixXi &FE,
                          const Eigen::MatrixXi &EV, const Eigen::MatrixX3d &B1,
                          const std::vector<int> &face_path,
                          const double &initial_angle,
                          std::vector<int> &edgeIds,
                          std::vector<int> &edgeSigns) {
  double transported_angle = initial_angle;

  for (int i = 1; i < face_path.size(); i++) {
    int prev_fid = face_path[i - 1];
    int post_fid = face_path[i];

    Eigen::Vector3d prevNormal = face_normal(V, F, prev_fid);
    Eigen::Vector3d postNormal = face_normal(V, F, post_fid);

    int edge_id = -1, edge_sign;

    get_edge_from_faces(prev_fid, post_fid, FE, EF, edge_id, edge_sign);

    Eigen::Vector3d edgeVec =
        (V.row(EV(edge_id, 1)).transpose() - V.row(EV(edge_id, 0)).transpose())
            .normalized();
    Eigen::Vector3d invertedEdgeVec =
        (V.row(EV(edge_id, 0)).transpose() - V.row(EV(edge_id, 1)).transpose())
            .normalized();

    double prev_angle =
        get_angle(B1.row(prev_fid).transpose(), edgeVec, prevNormal);
    double post_angle =
        get_angle(B1.row(post_fid).transpose(), edgeVec, postNormal);

    double edge_jump = -prev_angle + post_angle;

    double inverted_prev_angle =
        get_angle(B1.row(prev_fid).transpose(), invertedEdgeVec, prevNormal);
    double inverted_post_angle =
        get_angle(B1.row(post_fid).transpose(), invertedEdgeVec, postNormal);

    double inverted_edge_jump = -inverted_prev_angle + inverted_post_angle;

    edge_jump = std::abs(edge_jump) > std::abs(inverted_edge_jump) ? inverted_edge_jump : edge_jump;

    transported_angle = transported_angle - prev_angle + post_angle;

    edgeIds.emplace_back(edge_id);
    edgeSigns.emplace_back(edge_sign);
  }

  return transported_angle;
};
}  // namespace rp
