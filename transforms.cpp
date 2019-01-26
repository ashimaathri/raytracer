#include "transforms.h"

void Transform::rightmultiply(const glm::dmat4 &M, std::stack<glm::dmat4> &transforms) 
{
  glm::dmat4 &T = transforms.top();
  T = T * M;
}

glm::dmat4 Transform::rotate(const double degrees, const glm::dvec3& axis)
{
  const double angleInRadians = glm::radians(degrees);
  const double cosine = glm::cos(angleInRadians);
  const double sine = glm::sin(angleInRadians);

  const glm::dmat3 identity = glm::dmat3(1.0);

  const double x = axis.x;
  const double y = axis.y;
  const double z = axis.z;
  const glm::dmat3 dual = glm::transpose(glm::dmat3(0, -z, y, z, 0, -x, -y, x, 0));

  const glm::dmat3 rotationMatrix = cosine * identity + (1 - cosine) * outerProduct(axis, axis) + sine * dual;

  glm::dmat4 homogenizedRotation;
  homogenizedRotation[0] = glm::dvec4(rotationMatrix[0], 0);
  homogenizedRotation[1] = glm::dvec4(rotationMatrix[1], 0);
  homogenizedRotation[2] = glm::dvec4(rotationMatrix[2], 0);
  homogenizedRotation[3] = glm::dvec4(0, 0, 0, 1);

	return homogenizedRotation;
}

glm::dmat4 Transform::scale(const double &sx, const double &sy, const double &sz) 
{
  return glm::dmat4(sx, 0, 0, 0,
              0, sy, 0, 0,
              0, 0, sz, 0,
              0, 0, 0, 1);
}

glm::dmat4 Transform::translate(const double &tx, const double &ty, const double &tz) 
{
  return glm::dmat4(1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              tx, ty, tz, 1);
}
