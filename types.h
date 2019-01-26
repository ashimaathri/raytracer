#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stack>
#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <FreeImage.h>

#include "geometry.h"
#include "light.h"
#include "transforms.h"

#define MAX_LIGHTS 10
#define EPSILON 1e-4
#define MAX_SHAPES 100000

#ifndef _TYPES_H_
#define _TYPES_H_

namespace Type {
  class Ray;
  class AABB;
  class IntersectionInfo;
  class Shape;
  class Cube;
  class Sphere;
  class Triangle;
  class LineSegment;

  struct Material {
    glm::dvec3 ambient;
    glm::dvec3 diffuse;
    glm::dvec3 emission;
    glm::dvec3 specular;
    double shininess;
  };

  class Ray {
    public:
      glm::dvec3 origin;
      glm::dvec3 direction;

      Ray();
      Ray(glm::dvec3 o, glm::dvec3 d);

      Ray transform(glm::dmat4 transform);
      void print();

      glm::dvec3 getPoint(double distance);
  };

  class IntersectionInfo {
    public:
      double distance;
      glm::dvec3 point;
      glm::dvec3 displacedPoint;
      glm::dvec3 normal;
      Shape* shape;

      IntersectionInfo(double distance, glm::dvec3 point, glm::dvec3 normal, Shape* shape);
      IntersectionInfo();
  };

  class Shape {
    public:
      Material material;

      virtual bool does_intersect(Type::Ray ray, Type::IntersectionInfo &hit) = 0;
      virtual void print() = 0;
      virtual glm::dvec3 computeColorAtPoint(
          glm::dvec3 &eye,
          Type::IntersectionInfo &hit,
          glm::dvec3 &directionToLight,
          glm::dvec3 &lightColor);
  };

  class AABB : public Shape {
    public:
      double xmin;
      double ymin;
      double zmin;
      double xmax;
      double ymax;
      double zmax;
      Shape* shape;
      AABB* children[2];

      AABB() {};
      AABB(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, Shape* shape);
      Cube toCube();
      bool does_intersect(Type::Ray ray, Type::IntersectionInfo &hit);
      void print();
  };

  class LineSegment : public Shape {
    Ray ray;

    public:
      LineSegment();
      LineSegment(glm::dvec3 origin, glm::dvec3 direction);
      void print();
      glm::dvec3 computeColorAtPoint(
          glm::dvec3 &eye,
          Type::IntersectionInfo &hit,
          glm::dvec3 &directionToLight,
          glm::dvec3 &lightColor);
      bool does_intersect(Type::Ray ray, Type::IntersectionInfo &hit);
  };

  class Cube {
    public:
      static const int numCorners = 8;

      Cube(glm::dvec3 corners[numCorners]);
      Cube transform(glm::dmat4 transform);
      AABB* getAABB(Sphere* sphere);
      LineSegment** getEdges();

    private:
      glm::dvec3 corners[numCorners];
      AABB aabb;
  };

  class Sphere : public Shape {
    private:
      glm::dvec3 center;
      double radius;
      glm::dmat4 transform;
      glm::dmat4 inverseTransform;
      glm::dmat4 transposedInverseTransform;


    public:
      Sphere(double x, double y, double z, double r, glm::dvec3 ambient,
          glm::dvec3 diffuse, glm::dvec3 emission, glm::dvec3 specular,
          double shininess, glm::dmat4 transformationMatrix); 
      bool does_intersect(Type::Ray ray, Type::IntersectionInfo &hit);
      void print();
      AABB* getAABB();
  };

  class Triangle : public Shape {
    private:
      glm::dvec3* v1;
      glm::dvec3* v2;
      glm::dvec3* v3;
      glm::dvec3 normal;
      glm::dvec3 transformedNormal;
      glm::dmat4 transform;
      glm::dmat4 inverseTransform;
      glm::dmat4 transposedInverseTransform;

      glm::dvec3 computeNormals();
      bool barycentric_test(Ray ray, double distance);

    public:
      Triangle(glm::dvec3* vertex1, glm::dvec3* vertex2, glm::dvec3* vertex3,
          glm::dvec3 ambient, glm::dvec3 diffuse, glm::dvec3 emission,
          glm::dvec3 specular, double shininess, glm::dmat4 transformationMatrix);
      bool does_intersect(Type::Ray ray, Type::IntersectionInfo &hit);
      void print();
      AABB* getAABB();
  };
}
#endif
