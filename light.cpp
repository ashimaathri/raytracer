#include "light.h"

Light::DirectionalLight::DirectionalLight(glm::dvec3 direction, glm::dvec3 color) {
  this->direction = glm::normalize(direction);
  this->color = color;
}

glm::dvec3 Light::DirectionalLight::getColorAt(double distance) {
  return this->color;
}

glm::dvec3 Light::DirectionalLight::getDirectionFrom(glm::dvec3 point) {
  return this->direction;
}

double Light::DirectionalLight::getDistanceFrom(glm::dvec3 point) {
  return std::numeric_limits<double>::infinity();
}

Light::PointLight::PointLight(glm::dvec3 position, glm::dvec3 color, glm::dvec3 attenuation) {
  this->position = position;
  this->color = color;
  this->attenuation = attenuation;
}

glm::dvec3 Light::PointLight::getColorAt(double distance) {
  double denominator = glm::dot(this->attenuation, glm::dvec3(1, distance, std::pow(distance, 2)));
  return (1 /  denominator) * this->color;
}

glm::dvec3 Light::PointLight::getDirectionFrom(glm::dvec3 point) {
  return glm::normalize(this->position - point);
}

double Light::PointLight::getDistanceFrom(glm::dvec3 point) {
  return glm::length(this->position - point);
}
