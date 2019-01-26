#include <glm/glm.hpp>

#ifndef _LIGHTS_H_
#define _LIGHTS_H_

namespace Light {
  class Light;
  class DirectionalLight;
  class PointLight;

  class Light {
    public:
      virtual glm::dvec3 getColorAt(double distance) = 0;
      virtual glm::dvec3 getDirectionFrom(glm::dvec3 point) = 0;
      virtual double getDistanceFrom(glm::dvec3 point) = 0;
  };

  class DirectionalLight : public Light {
    private:
      glm::dvec3 color;
      glm::dvec3 direction;

    public:
      DirectionalLight() {};
      DirectionalLight(glm::dvec3 direction, glm::dvec3 color);
      glm::dvec3 getColorAt(double distance);
      glm::dvec3 getDirectionFrom(glm::dvec3 point);
      double getDistanceFrom(glm::dvec3 point);
  };

  class PointLight : public Light {
    private:
      glm::dvec3 color;
      glm::dvec3 position;
      glm::dvec3 attenuation;

    public:
      PointLight() {};
      PointLight(glm::dvec3 position, glm::dvec3 color, glm::dvec3 attenuation);
      glm::dvec3 getColorAt(double distance);
      glm::dvec3 getDirectionFrom(glm::dvec3 point);
      double getDistanceFrom(glm::dvec3 point);
  };
};

#endif
