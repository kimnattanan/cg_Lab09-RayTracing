#include "core/stats.h"
#include "core/vec3.h"
#include "student.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

// TASK 2a: Moller-Trumbore (DEF-10.01): does r hit the triangle at t > EPS? Write t, u, v.
bool ray_triangle(const Ray& r, const vec3& v0, const vec3& v1, const vec3& v2,
                  float& t, float& u, float& v) {
  vec3 e1 = v1-v0, e2 = v2-v0;
  vec3 h = cross(r.d,e2);
  float a = dot(e1,h);
  if(abs(a)<1e-8f) return false;
  vec3 s = r.o-v0, q = cross(s,e1);
  u = dot(s,h)/a;
  v = dot(r.d,q)/a;
  t = dot(e2,q)/a;
  return t>EPS && u>=0 && v>=0 && u+v<=1;
}

// TASK 2b: the closest triangle of m below t_max, testing every one (count_triangle_test()).
// Fill hit.t, hit.tri and hit.n.
bool intersect_mesh(const Ray& r, const Mesh& m, float t_max, Hit& hit) {
  bool ok = false;
  for(size_t i=0; i<m.tris.size(); ++i){
    count_triangle_test();
    const auto &[i0,i1,i2] = m.tris[i];
    float t,u,v;
    if(!ray_triangle(r, m.positions[i0], m.positions[i1], m.positions[i2], t, u, v) || t>t_max) continue;
    t_max = t;
    hit.t = t;
    hit.tri = i;
    hit.n = m.normals.empty() ?
      normalize(cross(m.positions[i1] - m.positions[i0], m.positions[i2] - m.positions[i0])) :
      normalize(m.normals[i0]*(1.f-u-v) + m.normals[i1]*u + m.normals[i2]*v);
    ok = true;
  }
  return ok;
}

// TASK 2c (optional): the slab test (EQ-10.03). inv_d is 1 / r.d; write the entry distance.
bool ray_aabb(const Ray& r, const vec3& inv_d, const AABB& box, float t_max, float& t_entry) {
  float t[2][3];
  for(int k=0; k<3; ++k){
    float t0 = (box.lo[k] - r.o[k]) * inv_d[k];
    float t1 = (box.hi[k] - r.o[k]) * inv_d[k];
    t[0][k] = std::min(t0,t1);
    t[1][k] = std::max(t0,t1);
  }
  t_entry = *std::max_element(t[0],t[0]+3);
  float t_exit = *std::min_element(t[1],t[1]+3);
  return t_entry < t_exit && t_exit > 0.f && t_entry < t_max;
}

// TASK 2d (optional): TASK 2b through m.bvh (count_node_visit() for every node).
bool intersect_bvh(const Ray& r, const Mesh& m, float t_max, Hit& hit) {
  if(m.bvh.nodes.empty()) return false;
  vec3 inv_d{1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z};
  int st[128];
  int sp = 0;
  st[sp++] = 0;
  bool ok = false;
  while(sp>0){
    int id = st[--sp];
    count_node_visit();
    const BVHNode &node = m.bvh.nodes[id];
    float t_entry;
    if(!ray_aabb(r, inv_d, node.box, t_max, t_entry)) continue;
    if(node.count){
      for(int i=0; i<node.count; ++i){
        int tri = m.bvh.order[node.first+i];
        count_triangle_test();
        const auto &[i0,i1,i2] = m.tris[tri];
        float t,u,v;
        if(!ray_triangle(r, m.positions[i0], m.positions[i1], m.positions[i2], t, u, v) || t>t_max) continue;
        t_max = t;
        hit.t = t;
        hit.tri = tri;
        hit.n = m.normals.empty() ?
          normalize(cross(m.positions[i1] - m.positions[i0], m.positions[i2] - m.positions[i0])) :
          normalize(m.normals[i0]*(1.f-u-v) + m.normals[i1]*u + m.normals[i2]*v);
        ok = true;
      }
    } else {
      int left = id+1, right = node.right;
      float t_l,t_r;
      bool hit_l = ray_aabb(r, inv_d, m.bvh.nodes[left].box, t_max, t_l);
      bool hit_r = ray_aabb(r, inv_d, m.bvh.nodes[right].box, t_max, t_r);
      if(hit_l && hit_r){
        if(t_l > t_r) std::swap(left,right);
        st[sp++] = left;
        st[sp++] = right;
      } else if(hit_l){
        st[sp++] = left;
      } else if(hit_r){
        st[sp++] = right;
      }
    }
  }
  return ok;
}
