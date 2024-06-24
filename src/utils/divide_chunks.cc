#include "divide_chunks.h"

namespace rp {
namespace divide_chunks {
namespace state {

class State {
 private:
  // どのsegmentにラベル付けされたか。デフォルト値は-1
  std::vector<int> face_id_label = {};
  // 今まで何個のfaceがラベル付けされたか
  int marked_faces_count = 0;
  // 今のセグメントの番号
  int current_segment_label = 0;
  // あるfaceにどのvertexが接続しているか
  std::vector<std::vector<int>> face_vertices = {};
  // あるvertexがどのfaceに接続しているか
  std::vector<std::vector<int>> vertex_faces = {};
  // あるfaceがどのfaceに接続しているか
  std::vector<std::vector<int>> face_faces = {};

 public:
  State(){};

  // mutations

  void initialize_state(int vertices_num, int faces_num) {
    face_id_label.resize(faces_num, -1);
    face_vertices.resize(faces_num, {});
    vertex_faces.resize(vertices_num, {});
    face_faces.resize(faces_num, {});
  };

  void register_face_and_vertex(int vid, int fid) {
    vertex_faces[vid].emplace_back(fid);
    face_vertices[fid].emplace_back(vid);
  }

  void register_face_faces(int fid1, int fid2) {
    face_faces[fid1].emplace_back(fid2);
    face_faces[fid2].emplace_back(fid1);
  }

  void compute_face_faces() {
    for (int i = 0; i < face_vertices.size(); i++) {
      auto current_fid = i;
      auto current_fid_vertices = face_vertices[i];

      for (int j = 0; j < current_fid_vertices.size(); j++) {
        auto vid = current_fid_vertices[j];
        auto faces_around_vids = vertex_faces[vid];

        for (int k = 0; k < faces_around_vids.size(); k++) {
          auto fid = faces_around_vids[k];

          // 同じものは加えない
          if (fid == current_fid) {
            continue;
          }

          // すでに入ってるものもここで飛ばす
          if (is_face_faces_registered(fid, current_fid)) {
            continue;
          }

          auto fid_vertices = face_vertices[fid];

          // fid_verticesとcurrent_fid_verticesの内2つが一致したらface_facesに追加する
          int count = 0;
          for (int l = 0; l < current_fid_vertices.size(); l++) {
            for (int m = 0; m < fid_vertices.size(); m++) {
              if (current_fid_vertices[l] == fid_vertices[m]) {
                count++;
              }
            }
          }

          if (count == 2) {
            register_face_faces(fid, current_fid);
          }
        }
      }
    }
  }

  void set_label(int fid) {
    if (face_id_label[fid] == -1) {
      face_id_label[fid] = current_segment_label;
      marked_faces_count++;
    }
  }

  void go_to_next_label() { current_segment_label++; }

  // getters

  bool is_face_faces_registered(int fid1, int fid2) {
    auto faces = face_faces[fid1];

    for (auto fid : faces) {
      if (fid == fid2) {
        return true;
      }
    }

    return false;
  }

  int get_new_segment_fid() {
    for (int i = 0; i < face_id_label.size(); i++) {
      if (face_id_label[i] == -1) {
        return i;
      }
    }

    return -1;
  };

  bool is_marked(int fid) { return face_id_label[fid] != -1; }

  bool is_labeling_finished() {
    return marked_faces_count >= face_id_label.size();
  }

  std::vector<int> get_adjacent_fids(int fid) { return face_faces[fid]; }

  int get_current_segment_label() { return current_segment_label; }

  std::vector<int> get_face_id_label() { return face_id_label; }

  int get_face_id_label_by_fid(int fid) { return face_id_label[fid]; }

  std::vector<std::vector<int>> get_face_faces() { return face_faces; }

  std::vector<std::vector<int>> get_vertex_faces() { return vertex_faces; }

  std::vector<int> get_vertex_face(int vid) { return vertex_faces[vid]; }

  std::vector<std::vector<int>> get_face_vertices() { return face_vertices; }
};

};  // namespace state

inline void initialize_state(Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
                             state::State &state) {
  state.initialize_state(V.rows(), F.rows());

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      state.register_face_and_vertex(F(i, j), i);
    }
  }

  state.compute_face_faces();
}

void recursive_labeling(std::vector<int> &fids, state::State &state) {
  // std::stringがあると末尾最適化が効かないことに注意
  // https://cpplover.blogspot.com/2018/04/c.html

  // なお、即時関数でのstd::coutは許容されるらしい。
  // [] { std::cout << "0000000";}() << std::endl;

  if (fids.size() == 0) return;

  std::vector<int> next_fids = {};

  for (int i = 0; i < fids.size(); i++) {
    auto fid = fids[i];

    if (state.is_marked(fids[i])) continue;

    state.set_label(fid);

    auto adjacent_fids = state.get_adjacent_fids(fid);

    for (int j = 0; j < adjacent_fids.size(); j++) {
      if (state.is_marked(adjacent_fids[j])) continue;

      next_fids.emplace_back(adjacent_fids[j]);
    }
  }

  recursive_labeling(next_fids, state);
}

inline void label_face_ids(state::State &state) {
  while (!state.is_labeling_finished()) {
    int start_id = state.get_new_segment_fid();

    std::cout << "recursive_labeling" << std::endl;

    std::vector<int> fids = {start_id};

    recursive_labeling(fids, state);

    state.go_to_next_label();
  }
}

inline void get_initial_border_point(
    Eigen::MatrixX3i &F, std::vector<std::vector<int>> &face_faces,
    std::vector<std::vector<int>> &vertex_faces,
    std::vector<std::vector<int>> &face_vertices,
    std::vector<std::vector<int>> &v_connection, int &face_id, int &vertex_id,
    int &next_vertex_id) {
  for (int i = 0; i < face_faces.size(); i++) {
    if (face_faces[i].size() != 2) {
      continue;
    }

    auto face = i;

    int inside_vid = -1;

    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        auto adjacent_face_1 = face_faces[i][0];
        auto adjacent_face_2 = face_faces[i][1];

        if (face_vertices[adjacent_face_1][j] ==
            face_vertices[adjacent_face_2][k]) {
          inside_vid = face_vertices[adjacent_face_1][j];
          break;
        }
      }
    }

    int border_vid = -1;
    int next_border_vid = -1;

    for (int j = 0; j < 3; j++) {
      auto vid = F(face, j);
      if (vid == inside_vid) {
        border_vid = F(face, (j + 1) % 3);
        next_border_vid = F(face, (j + 2) % 3);
      }
    }

    if (border_vid == -1) assert(false && "border vid is -1");

    auto connection = v_connection[border_vid];

    // ここに含まれるコネクションの内一つでもこのsegmentに含まれていないものがあればOK
    for (int j = 0; j < connection.size(); j++) {
      auto overborder_vid = connection[j];

      if (vertex_faces[overborder_vid].size() == 0) {
        face_id = face;
        vertex_id = border_vid;
        next_vertex_id = next_border_vid;
        return;
      }
    }
  }
};

inline void get_next_vertex_id(Eigen::MatrixX3i &F,
                               std::vector<std::vector<int>> &face_faces,
                               std::vector<std::vector<int>> &vertex_faces,
                               int &next_vid, int &next_fid) {
  int vertex_id = next_vid;
  // vertex_idと、その前のvertex_idの間にあるface
  int face_id = next_fid;
  auto candidate_faces = vertex_faces[vertex_id];

  // 今のface_idが2辺を境界に持つようなやつで、かつ今のvertex_idに接している三角形がこのface_idのみだったら、その三角形の次の点をnext_vidとする
  if (vertex_faces[vertex_id].size() == 1) {
    for (int i = 0; i < 3; i++) {
      if (F(face_id, i) == vertex_id) {
        next_vid = F(face_id, (i + 1) % 3);
        return;
      }
    }
  }

  for (int i = 0; i < candidate_faces.size(); i++) {
    auto candidate_fid = candidate_faces[i];
    auto adjacent_fids = face_faces[candidate_fid];

    // 周囲を全て囲まれているfaceは次のfaceとしては採用しない
    if (adjacent_fids.size() == 3) {
      continue;
    }

    // このfaceが、vertex_idを切り込みに含まない場合もダメ
    if (adjacent_fids.size() == 2) {
      auto adjacent_fid1 = adjacent_fids[0];
      auto adjacent_fid2 = adjacent_fids[1];
      int common_vid = -1;

      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          if (F(adjacent_fid1, j) == F(adjacent_fid2, k)) {
            common_vid = F(adjacent_fid1, j);
            break;
          }
        }
      }

      if (common_vid == vertex_id) {
        continue;
      }
    }

    for (int j = 0; j < 3; j++) {
      auto vid = F(candidate_fid, j);
      if (vid == vertex_id) {
        // 次のvidは反時計回りに回したもの
        next_vid = F(candidate_fid, (j + 1) % 3);
        next_fid = candidate_fid;
        return;
      }
    }
  }

  assert(false && "no border further from this point");
}

inline std::vector<std::vector<std::vector<int>>> get_connection(
    int segment_num, std::vector<std::vector<int>> &borders,
    std::vector<std::vector<int>> &v_connection, state::State &state) {
  int all_segments = state.get_current_segment_label();

  std::vector<std::vector<std::vector<int>>> result;
  result.resize(segment_num, {});

  for (int sid = 0; sid < segment_num; sid++) {
    // segment idごとに抽出

    result[sid].resize(borders.size(), {});

    for (int j = 0; j < borders.size(); j++) {
      auto border = borders[j];

      // -1で初期化
      result[sid][j].resize(border.size(), -1);
    }
  }

  for (int i = 0; i < borders.size(); i++) {
    auto border = borders[i];

    for (int j = 0; j < 0; j++) {
      auto border_vid = border[j];
      auto overborder_vids = v_connection[border_vid];

      for (int k = 0; k < overborder_vids.size(); k++) {
        auto overborder_vid = overborder_vids[k];

        // overborder_vidがどのsegmentに属するかを判定する
        auto fids = state.get_vertex_face(overborder_vid);

        assert(fids.size() > 0 && "size of fid must be over than 0");

        auto segment_id = state.get_face_id_label_by_fid(fids[0]);

        assert(segment_id != -1 && "segment id must not be -1");

        result[segment_id][i][j] = overborder_vid;
      }
    }
  }

  return result;
}

inline void divide_faces_and_vertices(
    Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
    // 特定のvidに対してどのvidが同じ点として対応しているかを示す
    std::vector<std::vector<int>> &v_connection,
    std::vector<Eigen::MatrixX3i> &result_faces,
    std::vector<Eigen::MatrixX3d> &result_vertices,
    std::vector<std::vector<std::vector<int>>> &result_borders,
    std::vector<std::vector<std::vector<std::vector<int>>>> &result_connection,
    std::vector<std::vector<int>> &initial_boundary_loops,
    state::State &state) {
  int segment_num = state.get_current_segment_label();

  result_faces.resize(segment_num, Eigen::MatrixX3i(0, 3));
  result_vertices.resize(segment_num, Eigen::MatrixX3d(V.rows(), 3));
  result_borders.resize(segment_num, {});
  result_connection.resize(segment_num, {});

  std::vector<int> segments_nums(segment_num, 0);

  auto face_id_label = state.get_face_id_label();

  std::vector<bool> is_in_boundary_loop(V.size(), false);

  for (int j = 0; j < initial_boundary_loops.size(); j++) {
    for (int k = 0; k < initial_boundary_loops[j].size(); k++) {
      int vid = initial_boundary_loops[j][k];
      is_in_boundary_loop[vid] = true;
    }
  }

  // 一回めの操作では領域確保を一気にやるための準備
  for (int i = 0; i < face_id_label.size(); i++) {
    auto segment_id = face_id_label[i];

    if (segment_id != -1) {
      segments_nums[segment_id]++;
    }
  }

  // 2回目では数に応じて配列をresizeする
  for (int i = 0; i < segment_num; i++) {
    auto &segment_faces = result_faces[i];
    segment_faces.resize(segments_nums[i], 3);

    result_vertices[i] = V;
  }

  std::vector<int> current_max_nums(segment_num, 0);

  // 3回目では代入
  for (int i = 0; i < face_id_label.size(); i++) {
    auto segment_id = face_id_label[i];

    assert(segment_id != -1 && "segment_id should not be -1");

    if (segment_id != -1) {
      result_faces[segment_id](current_max_nums[segment_id], 0) = F(i, 0);
      result_faces[segment_id](current_max_nums[segment_id], 1) = F(i, 1);
      result_faces[segment_id](current_max_nums[segment_id], 2) = F(i, 2);

      current_max_nums[segment_id]++;
    }
  }

  std::cout << segment_num << std::endl;

  for (int i = 0; i < segment_num; i++) {
    // 3. result_bordersとresult_connectionsを置き換える
    std::vector<std::vector<int>> loop;
    igl::boundary_loop(result_faces[i], loop);

    assert(loop.size() <= 2);

    // loop sizeが2なら、loop[0]が底面となるようにする
    if (loop.size() == 2) {
      bool no_invert_flag = false;

      // 底面を含むことの判定にinitial_boundary_loops（色々な処理を回す前のboundary_loop）を使う
      for (int j = 0; j < loop[0].size(); j++) {
        if (is_in_boundary_loop[loop[0][j]]) {
          no_invert_flag = true;
          break;
        }
      }

      if (!no_invert_flag) {
        auto tmp0 = loop[0];
        auto tmp1 = loop[1];

        loop[0] = tmp1;
        loop[1] = tmp0;
      }
    }

    result_borders[i] = loop;

    result_connection[i] =
        get_connection(segment_num, result_borders[i], v_connection, state);
  }
}

void remove_unreferenceds(
    std::vector<Eigen::MatrixX3i> &result_faces,
    std::vector<Eigen::MatrixX3d> &result_vertices,
    std::vector<std::vector<std::vector<int>>> &result_border,
    std::vector<std::vector<std::vector<std::vector<int>>>>
        &result_connection) {
  for (int sid = 0; sid < result_faces.size(); sid++) {
    auto F = result_faces[sid];
    auto V = result_vertices[sid];

    Eigen::MatrixX3i NF;
    Eigen::MatrixX3d NV;
    Eigen::MatrixXi I;

    igl::remove_unreferenced(V, F, NV, NF, I);

    // borderをupdateする
    for (int i = 0; i < result_border[sid].size(); i++) {
      for (int j = 0; j < result_border[sid][i].size(); j++) {
        result_border[sid][i][j] = I(result_border[sid][i][j]);
      }
    }

    // result_connectionをupdateする
    for (int i = 0; i < result_connection.size(); i++) {
      // 各result_connectionの内当該segment_idの中身を置き換える
      for (int j = 0; j < result_connection[i][sid].size(); j++) {
        for (int k = 0; k < result_connection[i][sid][j].size(); k++) {
          if (result_connection[i][sid][j][k] != -1) {
            result_connection[i][sid][j][k] =
                I(result_connection[i][sid][j][k]);
          }
        }
      }
    }

    result_faces[sid] = NF;
    result_vertices[sid] = NV;
  }
}

// つながっていない辺で分離してconnectionを作る。
void process(
    Eigen::MatrixX3d &V, Eigen::MatrixX3i &F,
    std::vector<std::vector<int>> &v_connection,
    std::vector<Eigen::MatrixX3i> &result_faces,
    std::vector<Eigen::MatrixX3d> &result_vertices,
    std::vector<std::vector<std::vector<int>>> &result_border,
    std::vector<std::vector<std::vector<std::vector<int>>>> &result_connection,
    std::vector<std::vector<int>> &initial_boundary_loops) {
  auto state = state::State();

  std::cout << "initialize vertex faces" << std::endl;

  // 0. stateの初期化
  initialize_state(V, F, state);

  std::cout << "label face ids" << std::endl;

  // 1. 全てのfaceをラベリングして分割
  label_face_ids(state);

  auto num = state.get_current_segment_label();
  std::cout << "current segment lable: " << num << std::endl;

  std::cout << "divide faces and vertices" << std::endl;

  // 2. 全てのfaceとvertexを配列に分割
  // result_border
  divide_faces_and_vertices(V, F, v_connection, result_faces, result_vertices,
                            result_border, result_connection,
                            initial_boundary_loops, state);

  // 3. unnnecesary verticesを削除する
  remove_unreferenceds(result_faces, result_vertices, result_border,
                       result_connection);
};

};  // namespace divide_chunks
};  // namespace rp
