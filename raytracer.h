#include "types.h"
#include "kdtree.h"
#include "omp.h"

namespace Raytracer {
  class Camera;

  class Camera {
    public:
      glm::dvec3 eye;
      glm::dvec3 lookAt;
      glm::dvec3 up;
      double fovy;
      glm::dmat3 orthonormalBasis;

      Camera() {};
      Camera(glm::dvec3 eye, glm::dvec3 lookAt, glm::dvec3 up, double fovy);
      void print();

    private:
      void computeOrthonormalBasis();
  };

  class Scene {
    public:
      int width;
      int height;
      int maxDepth;
      std::string outputFilename;
      Camera camera;
      int numVertices;
      glm::dvec3** vertices;
      int numLights;
      Light::Light* lights[MAX_LIGHTS];
      int numShapes;
      Type::Shape* shapes[MAX_SHAPES];
      Type::AABB* aabbs[MAX_SHAPES];
      Type::AABB* aabbHierarchy;
      std::stack<glm::dvec3> ambient;
      std::stack<glm::dvec3> diffuse;
      std::stack<glm::dvec3> emission;
      std::stack<glm::dvec3> specular;
      std::stack<double> shininess;
      std::stack<glm::dvec3> attenuation;
      std::stack<glm::dmat4> transforms;

      Scene(const char* filename);
      ~Scene();
      Type::IntersectionInfo closestIntersection(Type::Ray ray, Type::AABB* aabb);
      Type::IntersectionInfo anyIntersection(Type::Ray ray, Type::AABB* aabb,
          double maxDistance=std::numeric_limits<double>::infinity());
      void draw();
      void print();
      BYTE* render();
      bool isPointVisibleToLight(glm::dvec3 &hitPoint, Light::Light* light, glm::dvec3 &directionToLight, glm::dvec3 &color);
      glm::dvec3 findColor(Type::IntersectionInfo &hit, glm::dvec3 eye, int recurDepth);
      glm::dvec3 computeReflectionAtPoint(Type::IntersectionInfo &hit, glm::dvec3 &eye, int recurDepth);

      static void printProgress(int row, int column, int area, int& printed);
      static void setPixelColors(int row, int column, int width, BYTE* image, char color[3]);
      static char* colorToCharArray(glm::dvec3 color);
      static glm::dvec3 rayDirection(double row, double column, double half_width,
          double half_height, double tan_half_fovy, double tan_half_fovx,
          glm::dmat3 cameraBasis);
      static glm::dvec3 computeColorAtPoint(glm::dvec3 &eye,
          Type::IntersectionInfo &hit, glm::dvec3 &directionToLight,
          glm::dvec3 &color);
      static glm::dvec3 displacePointAlongNormal(glm::dvec3 &point, glm::dvec3 &normal);

    private:
      void processLine(std::string line);
      void parse(const char* filename);
  };
};
