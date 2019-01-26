#include <stack>
#include <glm/glm.hpp>

#ifndef _TRANSFORMS_H_
#define _TRANSFORMS_H_
namespace Transform {
  void rightmultiply(const glm::dmat4 &M, std::stack<glm::dmat4> &transforms);
  glm::dmat4 rotate(const double degrees, const glm::dvec3& axis);
  glm::dmat4 scale(const double &sx, const double &sy, const double &sz);
  glm::dmat4 translate(const double &tx, const double &ty, const double &tz);
}
#endif
