#include "core/scene.h"
#include "core/vec3.h"
#include "student.h"

#include <algorithm>
#include <cmath>

// TASK 1a: the ray from the camera through the centre of pixel (x, y).
Ray camera_ray(int x, int y, int width, int height, float fov) {
  float m = tan(fov/2.f);
  return Ray{vec3{0, 0, 0}, normalize(vec3{
    (2.f * (x+0.5f) / width - 1.f) * m * width / height,
    -(2.f * (y+0.5f) / height - 1.f) * m,
    -1.f
  })};
}

// TASK 1b: the smallest t > EPS where r meets the sphere (EQ-10.02).
bool ray_sphere(const Ray& r, const Sphere& s, float& t) {
  vec3 oc = r.o - s.center;
  float a = dot(r.d,r.d);
  float b = 2.f * dot(r.d,oc);
  float c = dot(oc,oc) - s.radius * s.radius;
  float disc = b*b - 4.f*a*c;
  if(disc < 0) return false;
  
  float rt = std::sqrt(disc);
  float t1 = (-b-rt)/(2.f*a);
  float t2 = (-b+rt)/(2.f*a);
  if(t1 > EPS) return t=t1, true;
  if(t2 > EPS) return t=t2, true;
  return false;
}

// TASK 1c: d reflected about the unit normal n (EQ-10.04).
vec3 reflect(const vec3& d, const vec3& n) {
  return d - n * (2.f * dot(d,n));
}

// TASK 1c: d refracted through a surface whose normal n points out of the object.
// false on total internal reflection.
bool refract(const vec3& d, const vec3& n, float ior, vec3& out) {
  float e1 = 1.f, e2 = ior;
  float cosi = std::clamp(dot(d,n), -1.f, 1.f);
  vec3 n2 = n;
  if(cosi < 0.f) cosi = -cosi;
  else n2 = -n, std::swap(e1,e2);

  float eta = e1/e2;
  float k = 1.f - eta * eta * (1.f - cosi * cosi);
  if(k < 0.f) return false;
  out = d * eta + n2 * (eta * cosi - std::sqrt(k));
  return true;
}

static vec3 calc(vec3 p,vec3 n,vec3 dir){
  return p + n * (dot(dir,n)<0.f ? -EPS : EPS);
}

// TASK 1d: the colour seen along r: Phong, shadows, reflection, refraction (Whitted).
vec3 cast_ray(const Ray& r, const Scene& scene, int depth) {
  if(depth > MAX_DEPTH) return scene.background;
  Hit hit;
  if(!scene_intersect(r,scene,hit)) return scene.background;

  const Material &mat = *hit.material;
  if(STAGE < 2) return mat.diffuse;

  float diff = 0.f, spec = 0.f;
  for(auto &light : scene.lights){
    vec3 l = normalize(light.position - hit.p);
    if(STAGE >= 4) {
      vec3 o = calc(hit.p, hit.n, l);
      Hit tmp;
      if(scene_intersect(Ray{o,l},scene,tmp) && tmp.t < length(light.position - o)) continue;
    }
    diff += light.intensity * std::max(0.f, dot(l,hit.n));
    if(STAGE >= 3) spec += light.intensity * std::pow(std::max(0.f, -dot(reflect(-l,hit.n), r.d)), mat.specular_exponent);
  }

  vec3 col = mat.diffuse * diff * mat.albedo[0];
  if(STAGE >= 3) col += vec3{1,1,1} * spec * mat.albedo[1];
  if(STAGE >= 5){
    vec3 r2 = normalize(reflect(r.d,hit.n));
    col += cast_ray(Ray{calc(hit.p, hit.n, r2), r2}, scene, depth+1) * mat.albedo[2];
  }
  if(STAGE >= 6){
    vec3 r2;
    if(refract(r.d,hit.n,mat.ior,r2)){
      r2 = normalize(r2);
      col += cast_ray(Ray{calc(hit.p, hit.n, r2), r2}, scene, depth+1) * mat.albedo[3];
    }
  }
  return col;
}
