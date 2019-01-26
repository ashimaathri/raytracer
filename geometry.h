#include <glm/glm.hpp>

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

namespace Geometry {
  glm::dvec3 dehomogenize(glm::dvec4 point);
  glm::dvec3 transformPoint(glm::dvec3 point, glm::dmat4 transform);
  glm::dvec3 transformDirection(glm::dvec3 direction, glm::dmat4 transform);
}

#endif
