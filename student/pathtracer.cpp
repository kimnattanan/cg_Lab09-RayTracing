#include "core/scene.h"
#include "core/vec3.h"
#include "student.h"

#include <algorithm>
#include <cmath>
#include <limits>

// TASK 3a: a unit direction around n with density cos(theta) / pi.
vec3 cosine_sample_hemisphere(const vec3& n, float u1, float u2) {
  float r = std::sqrt(u1), theta = 2*PI*u2;
  float x = r*std::cos(theta), y = r*std::sin(theta), z = std::sqrt(1-u1);
  vec3 t,b;
  make_basis(n,t,b);
  return normalize(x*t + y*b + z*n);
}

// TASK 3b: the light reaching hit straight from one random point on the lamp.
vec3 sample_light(const Hit& hit, const PtScene& scene, Rng& rng) {
  const PtLight &light = scene.light;
  vec3 q = light.corner + light.e1*rng.next() + light.e2*rng.next();
  vec3 o = hit.p + hit.n * EPS;
  
  // original
  vec3 to = q-o;
  
  // trap 1
  // vec3 to = q-hit.p;
  
  float dist = length(to);
  vec3 l = normalize(q-o);
  float cos_x = dot(l,hit.n);
  float cos_l = -dot(l,light.normal());
  if(cos_x <= 0.f || cos_l <= 0.f) return vec3{0,0,0};

  Hit tmp;
  if(pt_intersect(Ray{o,l}, scene, dist-EPS, tmp)) return vec3{0,0,0};

  const vec3 &rho = pt_material(scene,hit).albedo;
  float p = 1.f / light.area();
  return rho * (1.f/PI) * light.emission * (cos_x * cos_l / (dist*dist)) / p;
}

// TASK 3c: Russian roulette (EQ-10.08): false stops the path.
bool russian_roulette(vec3& throughput, int bounce, Rng& rng) {
  if(bounce < 3) return true;
  float q = std::clamp(1.f - luminance(throughput), 0.05f, 0.95f);
  if(rng.next() < q) return false;
  throughput /= (1.f-q);
  return true;
}

static const float INF = std::numeric_limits<float>::infinity();

// TASK 3d: the light arriving along r, following one path.
vec3 path_trace(const Ray& r, const PtScene& scene, Rng& rng, const PtOptions& opt) {
  vec3 throughput{1,1,1};
  vec3 res{0,0,0};
  Ray ray = r;
  for(int bounce=0; bounce<opt.max_depth; ++bounce){
    Hit hit;
    if(!pt_intersect(ray,scene,INF,hit)) break;
    if(dot(hit.n,ray.d) > 0) hit.n = -hit.n;

    const PtMaterial &mat = pt_material(scene,hit);

    if(luminance(mat.emission) > 0.f){
      // original
      if((bounce==0 || !opt.nee) && dot(ray.d,scene.light.normal()) < 0.f)
        res += throughput * mat.emission;
      
      // trap 2
      // if((bounce==0 || !opt.nee))
      //   res += throughput * mat.emission;
      
      break;
    }
    if(opt.nee) res += throughput * sample_light(hit,scene,rng);
    vec3 dir = cosine_sample_hemisphere(hit.n, rng.next(), rng.next());
    throughput *= mat.albedo;
    ray = Ray{hit.p + hit.n*EPS, dir};
    
    if(opt.rr && !russian_roulette(throughput,bounce,rng)) break;
  }
  return res;
}
