#include "geometry.h"

glm::dvec3 Geometry::dehomogenize(glm::dvec4 point) {
  return glm::dvec3(point) / point.w;
}

glm::dvec3 Geometry::transformPoint(glm::dvec3 point, glm::dmat4 transform) {
  return dehomogenize(transform * glm::dvec4(point, 1));
}

glm::dvec3 Geometry::transformDirection(glm::dvec3 direction, glm::dmat4 transform) {
  return glm::dvec3(transform * glm::dvec4(direction, 0));
}
