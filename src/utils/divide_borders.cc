#include "divide_borders.h"

namespace rp {
namespace divide_borders {
namespace state {

struct FaceVertex {
  int face_id;
  int old_vertex_id;
  int new_vertex_id;

  FaceVertex(int _fid, int _ovid, int _nvid) {
    face_id = _fid;
    old_vertex_id = _ovid;
    new_vertex_id = _nvid;
  }
};

struct NewVertex {
  int vertex_id;
  int original_vertex_id;

  NewVertex(int _vid, int _ovid) {
    vertex_id = _vid;
    original_vertex_id = _ovid;
  }
};

struct VertexFaces {
  int vertex_id;
  // 当該vertexに接続するfaceの一覧。
  std::vector<int> face_ids;

  VertexFaces(int id) {
    vertex_id = id;
    face_ids = {};
  }

  VertexFaces() {}
};

struct VertexConnection {
  int vertex_id;
  // 当該vertexに接続するvertexの一覧。
  std::vector<int> vertex_ids;
  // すでに点を展開する作業が行われたかどうか
  bool processed_flag = false;

  VertexConnection(int id) {
    vertex_id = id;
    vertex_ids = {};
  }

  VertexConnection() {}
};

class State {
 private:
  std::vector<FaceVertex> updating_face_vertices_list = {};
  std::vector<NewVertex> new_vertices_list = {};
  std::vector<VertexFaces> vertex_faces = {};
  std::vector<VertexConnection> vertex_connections = {};
  int current_max_index = 0;

 public:
  // mutations
  void initialize_vertex_faces(int num) {
    vertex_faces.clear();

    vertex_faces.resize(num);

    for (int i = 0; i < num; i++) {
      vertex_faces[i] = state::VertexFaces(i);
    }
  }

  void add_face_to_vertex_faces(int vid, int fid) {
    vertex_faces[vid].face_ids.emplace_back(fid);
  }

  void initialize_vertex_connections(int vertex_nums) {
    vertex_connections.resize(vertex_nums);
    current_max_index = vertex_nums - 1;

    for (int i = 0; i < vertex_nums; i++) {
      vertex_connections[i] = VertexConnection(i);
    }
  }

  bool has_duplicate_vertex_connection(int vid, int new_vid) {
    // 重複があるか確認
    for (int i = 0; i < vertex_connections[vid].vertex_ids.size(); i++) {
      if (vertex_connections[vid].vertex_ids[i] == new_vid) {
        return true;
      }
    }

    return false;
  }

  void add_vertex_connection(int vid_1, int vid_2) {
    if (vid_1 != vid_2) {
      if (!has_duplicate_vertex_connection(vid_1, vid_2)) {
        vertex_connections[vid_1].vertex_ids.emplace_back(vid_2);
      }

      if (!has_duplicate_vertex_connection(vid_2, vid_1)) {
        vertex_connections[vid_2].vertex_ids.emplace_back(vid_1);
      }
    }
  }

  void toggle_vertex_connection_flag(int vid) {
    vertex_connections[vid].processed_flag = true;
  }

  void add_new_vertex(std::vector<int> &face_ids, int original_vertex_id) {
    current_max_index++;
    auto new_vertex_id = current_max_index;

    for (int i = 0; i < face_ids.size(); i++) {
      auto face_id = face_ids[i];

      auto fv = FaceVertex(face_id, original_vertex_id, new_vertex_id);
      updating_face_vertices_list.emplace_back(fv);

      auto nv = NewVertex(new_vertex_id, original_vertex_id);
      new_vertices_list.emplace_back(nv);
    }
  }

  // getters
  VertexConnection get_vertex_connection(int vid) {
    auto vc = vertex_connections[vid];
    return vc;
  }

  VertexFaces get_vertex_faces(int vid) {
    auto vf = vertex_faces[vid];
    return vf;
  }

  std::vector<FaceVertex> get_all_updating_face_vertices() {
    return updating_face_vertices_list;
  }

  std::vector<NewVertex> get_all_new_vertices() { return new_vertices_list; }
};
};  // namespace state

inline void initialize_vertex_faces(const Eigen::MatrixX3d &V,
                                    const Eigen::MatrixX3i &F,
                                    state::State &state) {
  int v_rows = V.rows();
  state.initialize_vertex_faces(v_rows);

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < F.cols(); j++) {
      auto vid = F(i, j);

      assert(vid < v_rows && "this vid is outside of vertex list");

      state.add_face_to_vertex_faces(vid, i);
    }
  }
}

inline void initialize_vertex_connections(
    const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F,
    const std::vector<Eigen::VectorXi> &positions, state::State &state) {
  state.initialize_vertex_connections(V.rows());

  for (int i = 0; i < positions.size(); i++) {
    auto path = positions[i];

    for (int j = 0; j < path.size(); j++) {
      if (j != 0) {
        state.add_vertex_connection(path(j), path(j - 1));
      }
    }
  }

  // boundaryもvertex_connectionsに登録する
  std::vector<std::vector<int>> boundary_loops;
  igl::boundary_loop(F, boundary_loops);

  for (int i = 0; i < boundary_loops.size(); i++) {
    auto path = boundary_loops[i];

    for (int j = 0; j < path.size(); j++) {
      state.add_vertex_connection(path[j],
                                  path[(j - 1 + path.size()) % path.size()]);
    }
  }
}

// kとlのindex以外のadjacent_vidsが当該faceに含まれるかどうか
inline bool contains_vertex(const int &current_face, const Eigen::MatrixX3i &F,
                            const std::vector<int> &adjacent_vids, const int &k,
                            const int &l) {
  for (int i = 0; i < 3; i++) {
    auto vid = F(current_face, i);

    for (int j = 0; j < adjacent_vids.size(); j++) {
      if (k == j || l == j) {
        continue;
      }

      if (adjacent_vids[j] == vid) {
        return true;
      }
    }
  }

  return false;
}

// center_vidと、current_faceの点の内center_vidでもfrontier_vidでもない点の両方を含む三角形を探す
inline int next_face(const int &current_face_id, const Eigen::MatrixX3i &F,
                     const int &center_vid, int &frontier_vid,
                     const std::vector<int> already_added_fids,
                     state::State &state) {
  auto fids = state.get_vertex_faces(center_vid).face_ids;

  int other_vid = -1;

  for (int i = 0; i < 3; i++) {
    auto current_vid = F(current_face_id, i);
    if (current_vid != center_vid && current_vid != frontier_vid) {
      other_vid = current_vid;
      break;
    }
  }

  assert(other_vid != -1 && "other vid not found in this face");

  for (int i = 0; i < fids.size(); i++) {
    auto fid = fids[i];

    if (fid == current_face_id) {
      continue;
    }

    if ((F(fid, 0) == other_vid || F(fid, 1) == other_vid ||
         F(fid, 2) == other_vid)) {
      bool continue_flag = false;

      for (int j = 0; j < already_added_fids.size(); j++) {
        if (already_added_fids[j] == fid) {
          continue_flag = true;
          break;
        }
      }

      if (continue_flag) {
        continue;
      } else {
        frontier_vid = other_vid;
        return fid;
      }
    }
  }

  assert(false && "next face is not found");
  return -1;
}

inline void execute_division(const Eigen::MatrixX3i &F,
                             const std::vector<Eigen::VectorXi> &positions,
                             state::State &state) {
  for (int i = 0; i < positions.size(); i++) {
    auto path = positions[i];

    for (int j = 0; j < path.size(); j++) {
      // path上の点でイテレーションを回す

      // path上での点p_jの隣接関係
      auto connection = state.get_vertex_connection(path(j));

      if (connection.processed_flag == false) {
        // path上で点p_jに隣接している点
        auto adjacent_vids = connection.vertex_ids;

        // 隣接点の内2つ（p_k, p_l）を選び、(p_k, p_j)と(p_l,
        // p_j)に囲まれた面は全て辺の番号を新たに作って書き換える。 なお、(p_k,
        // p_l)と(p_l,
        // p_j)の間に、adjacent_vidsに含まれる別の点が存在した場合には置き換えない。
        for (int k = 0; k < adjacent_vids.size(); k++) {
          for (int l = 0; l < k; l++) {
            auto p_j = path(j);
            auto p_k = adjacent_vids[k];
            auto p_l = adjacent_vids[l];
            const auto vertex_faces = state.get_vertex_faces(p_j);

            auto fids = vertex_faces.face_ids;
            // 辺(p_k, p_j)に接する2つの面のid
            std::vector<int> start_face_ids = {};

            for (int m = 0; m < fids.size(); m++) {
              auto fid = fids[m];

              if (F(fid, 0) == p_k || F(fid, 1) == p_k || F(fid, 2) == p_k) {
                start_face_ids.emplace_back(fid);
              }
            }

            assert(start_face_ids.size() > 0 && start_face_ids.size() < 3 &&
                   "start face id seems to be invalid");

            // (p_k, p_j)と(p_l, p_j)に囲まれた面のid
            std::vector<int> between_face_ids = {};

            for (auto start_face_id : start_face_ids) {
              // 辺(p_k, p_j)に接する面から開始し、隣接しかつp_jを含む面を辿っていく。

              // 今探索してるface
              int current_face = start_face_id;

              between_face_ids.clear();

              int frontier_vid = p_k;

              while (true) {
                if (F(current_face, 0) == p_l || F(current_face, 1) == p_l ||
                    F(current_face, 2) == p_l) {
                  // 今のfaceがp_lを含んでいる
                  between_face_ids.emplace_back(current_face);
                  break;
                } else if (contains_vertex(current_face, F, adjacent_vids, k,
                                           l)) {
                  // 今のfaceがp_lとp_k以外のadjacent_vidsを含んでいる場合
                  // clearする
                  between_face_ids.clear();
                  between_face_ids.resize(0);
                  break;
                } else {
                  // 次に行く
                  between_face_ids.emplace_back(current_face);

                  // current_faceと、p_jとfrontier_vidを共有する辺で接する面を次の面とする
                  current_face = next_face(current_face, F, p_j, frontier_vid,
                                           between_face_ids, state);
                }
              };

              if (between_face_ids.size() > 0) {
                // もしbetween_face_idsがOKならここで探索は終了
                break;
              }
            }

            state.add_new_vertex(between_face_ids, p_j);
          }
        }

        state.toggle_vertex_connection_flag(path(j));
      }
    }
  }
};

inline void update_mesh(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                        state::State &state,
                        std::vector<std::vector<int>> &v_connection,
                        std::vector<std::array<int, 2>> &new_vertices_pairs) {
  auto face_vertices = state.get_all_updating_face_vertices();

  // faceのvertex idを書き換え
  for (int i = 0; i < face_vertices.size(); i++) {
    auto fid = face_vertices[i].face_id;

    for (int j = 0; j < 3; j++) {
      if (F(fid, j) == face_vertices[i].old_vertex_id) {
        F(fid, j) = face_vertices[i].new_vertex_id;
      }
    }
  }

  auto new_vertices = state.get_all_new_vertices();

  // vertexを追加
  V.conservativeResize(V.rows() + new_vertices.size(), 3);
  // v_connectionは初期化する
  v_connection.resize(V.rows(), {});
  for (int i = 0; i < new_vertices.size(); i++) {
    auto original_vid = new_vertices[i].original_vertex_id;
    auto new_vid = new_vertices[i].vertex_id;

    V(new_vid, 0) = V(original_vid, 0);
    V(new_vid, 1) = V(original_vid, 1);
    V(new_vid, 2) = V(original_vid, 2);

    v_connection[original_vid].emplace_back(new_vid);
    v_connection[new_vid].emplace_back(original_vid);

    std::array<int, 2> nvp;
    nvp[0] = original_vid;
    nvp[1] = new_vid;
    new_vertices_pairs.emplace_back(nvp);
  }
};

void process(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
    // 切り開く点のindices。ループの場合は最後の要素と最初の要素が一致するようにする。
    const std::vector<std::vector<int>> &positions) {
  std::vector<std::vector<int>> v_connection;
  std::vector<std::array<int, 2>> new_vertices;

  std::vector<Eigen::VectorXi> positions_v;
  positions_v.resize(positions.size(), Eigen::VectorXi(0));

  for (int i = 0; i < positions.size(); i++) {
    Eigen::VectorXi vec;
    vec.resize(positions[i].size());
    for (int j = 0; j < positions[i].size(); j++) {
      vec(j) = positions[i][j];
    }
    positions_v[i] = vec;
  }

  process(V, F, positions_v, v_connection, new_vertices);
}

void process(
    Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
    // 切り開く点のindices。ループの場合は最後の要素と最初の要素が一致するようにする。
    const std::vector<Eigen::VectorXi> &positions,
    std::vector<std::vector<int>> &v_connection) {
  std::vector<std::array<int, 2>> new_vertices;
  process(V, F, positions, v_connection, new_vertices);
};

void process(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
             // 切り開く点のindices
             const std::vector<Eigen::VectorXi> &positions,
             std::vector<std::vector<int>> &v_connection,
             std::vector<std::array<int, 2>> &new_vertices) {
  auto state = state::State();

  // 1. 全てのface/vertexに対して接続関係を明らかにする
  initialize_vertex_faces(V, F, state);

  initialize_vertex_connections(V, F, positions, state);

  // 2. 各position上の点について分割を実行する
  execute_division(F, positions, state);

  // 3. VとFを書き換え、v_connectionに今回分離されたvertexを記録する
  update_mesh(V, F, state, v_connection, new_vertices);
};

};  // namespace divide_borders
};  // namespace rp
