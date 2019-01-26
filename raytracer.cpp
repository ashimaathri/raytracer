#include "raytracer.h"

Type::IntersectionInfo Raytracer::Scene::anyIntersection(Type::Ray ray, Type::AABB* aabb, double maxDistance) {
  Type::IntersectionInfo hit, emptyHit;

  if(!aabb->does_intersect(ray, hit)) {
    return hit;
  }

  if(aabb->shape) {
    aabb->shape->does_intersect(ray, hit);
    return hit.distance < maxDistance ? hit : emptyHit;
  }

  Type::IntersectionInfo leftHit = this->anyIntersection(ray, aabb->children[0], maxDistance);
  Type::IntersectionInfo rightHit = this->anyIntersection(ray, aabb->children[1], maxDistance);

  if(leftHit.shape && rightHit.shape) {
    return leftHit.distance < rightHit.distance ? leftHit : rightHit;
  } 

  return leftHit.shape ? leftHit : rightHit;
}

Type::IntersectionInfo Raytracer::Scene::closestIntersection(Type::Ray ray, Type::AABB* aabb) {
  Type::IntersectionInfo hit;

  if(!aabb->does_intersect(ray, hit)) {
    return hit;
  }

  if(aabb->shape) {
    aabb->shape->does_intersect(ray, hit);
    return hit;
  }

  Type::IntersectionInfo leftHit = this->closestIntersection(ray, aabb->children[0]);
  Type::IntersectionInfo rightHit = this->closestIntersection(ray, aabb->children[1]);

  if(leftHit.shape && rightHit.shape) {
    return leftHit.distance < rightHit.distance ? leftHit : rightHit;
  } 

  return leftHit.shape ? leftHit : rightHit;
}

void printHierarchy(Type::AABB* aabb) {
  std::cout << "node\n";
  if(aabb->shape) {
    aabb->shape->print();
    return;
  }
  if(aabb->children[0]) {
    std::cout << "left\n";
    printHierarchy(aabb->children[0]);
    std::cout << "right\n";
    printHierarchy(aabb->children[1]);
  }
}

Raytracer::Scene::Scene(const char* filename) : width(160), height(120), maxDepth(5), outputFilename("raytrace.png"), \
            numVertices(0), numLights(0), numShapes(0) {
  glm::dvec3 gray = glm::dvec3(0.2, 0.2, 0.2);
  glm::dvec3 black = glm::dvec3(0, 0, 0);

  this->ambient.push(gray);
  this->diffuse.push(black);
  this->emission.push(black);
  this->specular.push(black);
  this->shininess.push(0);
  this->attenuation.push(glm::dvec3(1, 0, 0));
  this->transforms.push(glm::dmat4(1.0f));
  this->parse(filename);
  this->aabbHierarchy = KDTree::buildKdTree(this->aabbs, 0, this->numShapes, 0);
  //printHierarchy(this->aabbHierarchy);
}

Raytracer::Scene::~Scene() {
  for(int i = 0; i < this->numLights; i++) {
    delete this->lights[i];
  }
  for(int i = 0; i < this->numShapes; i++) {
    delete this->shapes[i];
  }
}

void Raytracer::Scene::print() {
  std::cout << "Width: " << this->width << "\n";
  std::cout << "Height: " << this->height << "\n";
  std::cout << "MaxDepth: " << this->maxDepth << "\n";
  std::cout << "Output Filename: " << this->outputFilename << "\n";
  this->camera.print();
  for(int i = 0; i < this->numShapes; i++) {
    this->shapes[i]->print();
  }
  for(int i = 0; i < this->numVertices; i++) {
    std::cout << glm::to_string(*this->vertices[i]) << "\n";
  }
}

void Raytracer::Scene::draw() {
  FreeImage_Initialise();
  FIBITMAP *img = FreeImage_ConvertFromRawBits(this->render(),
                                               this->width,
                                               this->height,
                                               this->width * 3,
                                               24,
                                               0xFF0000,
                                               0x00FF00,
                                               0x0000FF,
                                               true);
  FreeImage_Save(FIF_PNG, img, this->outputFilename.c_str(), 0);
  FreeImage_DeInitialise();
}

void Raytracer::Scene::printProgress(int row, int column, int area, int& printed) {
  int percent_area_covered = (row * column * 100) / area;

  if(percent_area_covered % 10 == 0) {
    if(percent_area_covered > printed) {
      std::cout << percent_area_covered << " percent done \n";
      printed = percent_area_covered;
    }
  }
}

void Raytracer::Scene::setPixelColors(int row, int column, int width, BYTE* image, char color[3]) {
  static int numChannels = 3;

  for(int channel = 0; channel < numChannels; channel++) {
    int flatIndex = row * width * numChannels + column * numChannels + channel;
    image[flatIndex] = color[2 - channel];
  }
}

char* Raytracer::Scene::colorToCharArray(glm::dvec3 color) {
  char* byteArray = new char [3];

  glm::dvec3 clampedColor = glm::clamp(color, glm::dvec3(0), glm::dvec3(1));
  glm::dvec3 denormalizedColor = clampedColor * 255.0;

  byteArray[0] = denormalizedColor[0];
  byteArray[1] = denormalizedColor[1];
  byteArray[2] = denormalizedColor[2];

  return byteArray;
}

glm::dvec3 Raytracer::Scene::rayDirection(double row,
                       double column,
                       double half_width,
                       double half_height,
                       double tan_half_fovy,
                       double tan_half_fovx,
                       glm::dmat3 cameraBasis) {
  double alpha = tan_half_fovx * (column - half_width) / half_width;
  double beta = tan_half_fovy * (half_height - row) / half_height;
  return glm::normalize(alpha * cameraBasis[0] + beta * cameraBasis[1] - cameraBasis[2]);
}

BYTE* Raytracer::Scene::render() {
  const double half_width = this->width / 2;
  const double half_height = this->height / 2;
  const double aspect = (double)this->width / (double)this->height;
  const double tan_half_fovy = glm::tan(glm::radians(this->camera.fovy) / 2);
  const double tan_half_fovx = tan_half_fovy * aspect;
  const int area = this->width * this->height;
  int printed = 0;

  BYTE* image = new BYTE [this->height * this->width * 3];

#pragma omp parallel
#pragma omp for
  for(int row = 0; row < this->height; row++) {
    for(int column = 0; column < this->width; column++) {
      Raytracer::Scene::printProgress(row, column, area, printed);

      Type::Ray ray = Type::Ray(
          this->camera.eye,
          Raytracer::Scene::rayDirection(row + 0.5,
            column + 0.5,
            half_width,
            half_height,
            tan_half_fovy,
            tan_half_fovx,
            this->camera.orthonormalBasis));

      Type::IntersectionInfo hit = this->closestIntersection(ray, this->aabbHierarchy);

      Raytracer::Scene::setPixelColors(row,
          column,
          this->width,
          image,
          Raytracer::Scene::colorToCharArray(this->findColor(hit, this->camera.eye, 0)));
    }
  }

  return image;
}

bool Raytracer::Scene::isPointVisibleToLight(glm::dvec3 &hitPoint, Light::Light* light, glm::dvec3 &directionToLight, glm::dvec3 &lightColor) {
  double distanceToLight = light->getDistanceFrom(hitPoint);
  directionToLight = light->getDirectionFrom(hitPoint);
  lightColor = light->getColorAt(distanceToLight);

  Type::Ray shadowRay = Type::Ray(hitPoint, directionToLight);
  Type::IntersectionInfo shadowHit = this->anyIntersection(shadowRay, this->aabbHierarchy, distanceToLight);

  return !shadowHit.shape;
}

glm::dvec3 Raytracer::Scene::findColor(
    Type::IntersectionInfo &hit,
    glm::dvec3 eye,
    int recurDepth)
{
  if(!hit.shape) {
    return glm::dvec3(0, 0, 0);
  }
  glm::dvec3 lightColor;
  glm::dvec3 directionToLight;
  Type::Material material = hit.shape->material;
  glm::dvec3 computedColor = material.ambient + material.emission;

  for(int i = 0; i < this->numLights; i++) {
    if(this->isPointVisibleToLight(hit.displacedPoint, this->lights[i], directionToLight, lightColor)) {
      computedColor += hit.shape->computeColorAtPoint(eye, hit, directionToLight, lightColor);
    }
  }
  if(recurDepth == this->maxDepth) {
    return computedColor;
  }
  if(glm::all(glm::equal(material.specular, glm::dvec3(0)))) {
    return computedColor;
  }
  computedColor += this->computeReflectionAtPoint(hit, eye, recurDepth);

  return computedColor;
}

glm::dvec3 Raytracer::Scene::computeReflectionAtPoint(Type::IntersectionInfo &hit, glm::dvec3 &eye, int recurDepth) {
  glm::dvec3 directionToEye = glm::normalize(eye - hit.displacedPoint);
  double lDotN = glm::dot(directionToEye, hit.normal);
  glm::dvec3 mirrorDirection = glm::normalize(-directionToEye + 2.0 * lDotN * hit.normal);

  Type::Ray reflectedRay = Type::Ray(hit.displacedPoint, mirrorDirection);
  Type::IntersectionInfo reflectionHit = this->closestIntersection(reflectedRay, this->aabbHierarchy);
  glm::dvec3 specular = hit.shape->material.specular;
  return specular * this->findColor(reflectionHit, hit.displacedPoint, recurDepth + 1);
}

const bool isBlank(std::string line) {
  return line.find_first_not_of(" \t\r\n") == std::string::npos;
}

const bool isComment(std::string line) {
  return line[0] == '#';
}

bool readvals(std::stringstream &s, const int numvals, double* values) {
  for (int i = 0; i < numvals; i++) {
    s >> values[i];
    if (s.fail()) {
      std::cout << "Failed reading value " << i << " will skip\n";
      return false;
    }
  }
  return true; 
}

void Raytracer::Scene::parse(const char* filename) {
  std::ifstream in;
  in.open(filename);

  if (in.is_open()) {
    std::string line;
    getline (in, line);

    while (in) {
      this->processLine(line);
      getline (in, line); 
    }
  } else {
    std::cerr << "Unable to Open Input Scene File " << filename << "\n"; 
    exit(1); 
  }
}

void Raytracer::Scene::processLine(std::string line) {
  double values[10];

  if (isBlank(line) || isComment(line)) {
    return;
  }

  std::stringstream lineStream(line);

  std::string cmd;
  lineStream >> cmd;

  if (cmd == "size") {
    if(readvals(lineStream, 2, values)) {
      this->width = values[0];
      this->height = values[1];
    }
  }
  else if (cmd == "maxdepth") {
    if(readvals(lineStream, 1, values)) {
      this->maxDepth = values[0];
    }
  }
  else if (cmd == "output") {
    lineStream >> this->outputFilename;
    if (lineStream.fail()) {
      std::cout << "Failed reading value " << "output" << " will skip\n";
    }
  }
  else if (cmd == "camera") {
    if(readvals(lineStream, 10, values)) {
      glm::dvec3 eye = glm::dvec3(values[0], values[1], values[2]);
      glm::dvec3 lookAt = glm::dvec3(values[3], values[4], values[5]);
      glm::dvec3 up = glm::dvec3(values[6], values[7], values[8]);
      double fovy = values[9];
      this->camera = Raytracer::Camera(eye, lookAt, up, fovy);
    }
  }
  else if (cmd == "sphere") {
    if(readvals(lineStream, 4, values)) {
      Type::Sphere* sphere = new Type::Sphere(values[0], values[1], values[2],
          values[3], this->ambient.top(), this->diffuse.top(),
          this->emission.top(), this->specular.top(), this->shininess.top(),
          this->transforms.top());
      this->shapes[this->numShapes] = sphere;
      this->aabbs[this->numShapes] = sphere->getAABB();
      this->numShapes++;
      /*
      Type::LineSegment** aabbEdges = sphere->aabb.toCube().getEdges();
      for(int i = 0; i < 12; i++) {
        this->shapes[this->numShapes++] = aabbEdges[i];
      }
      */
    }
  }
  else if (cmd == "maxverts") {
    if(readvals(lineStream, 1, values)) {
      this->vertices = new glm::dvec3* [(int)values[0]];
    }
  }
  else if (cmd == "vertex") {
    if(readvals(lineStream, 3, values)) {
      this->vertices[this->numVertices++] = new glm::dvec3(values[0], values[1], values[2]);
    }
  }
  else if (cmd == "tri") {
    if(readvals(lineStream, 3, values)) {
      glm::dvec3* vertices[3];

      for(int i = 0; i < 3; i++) {
        vertices[i] = this->vertices[(int)values[i]];
      }

      Type::Triangle* triangle = new Type::Triangle(vertices[0], vertices[1],
          vertices[2], this->ambient.top(), this->diffuse.top(),
          this->emission.top(), this->specular.top(), this->shininess.top(),
          this->transforms.top());
      this->shapes[this->numShapes] = triangle;
      this->aabbs[this->numShapes] = triangle->getAABB();
      this->numShapes++;
      /*
      Type::LineSegment** aabbEdges = triangle->aabb.toCube().getEdges();
      for(int i = 0; i < 12; i++) {
        this->shapes[this->numShapes++] = aabbEdges[i];
      }
      */
    }
  }
  else if (cmd == "ambient") {
    if(readvals(lineStream, 3, values)) {
      this->ambient.push(glm::dvec3(values[0], values[1], values[2]));
    }
  }
  else if (cmd == "diffuse") {
    if(readvals(lineStream, 3, values)) {
      this->diffuse.push(glm::dvec3(values[0], values[1], values[2]));
    }
  }
  else if (cmd == "emission") {
    if(readvals(lineStream, 3, values)) {
      this->emission.push(glm::dvec3(values[0], values[1], values[2]));
    }
  }
  else if (cmd == "specular") {
    if(readvals(lineStream, 3, values)) {
      this->specular.push(glm::dvec3(values[0], values[1], values[2]));
    }
  }
  else if (cmd == "shininess") {
    if(readvals(lineStream, 1, values)) {
      this->shininess.push(values[0]);
    }
  }
  else if (cmd == "attenuation") {
    if(readvals(lineStream, 3, values)) {
      this->attenuation.push(glm::dvec3(values[0], values[1], values[2]));
    }
  }
  else if (cmd == "directional") {
    if(readvals(lineStream, 6, values)) {
      glm::dvec3 direction = glm::dvec3(values[0], values[1], values[2]);
      glm::dvec3 color = glm::dvec3(values[3], values[4], values[5]);

      this->lights[this->numLights++] = new Light::DirectionalLight(direction, color);
    }
  }
  else if (cmd == "point") {
    if(readvals(lineStream, 6, values)) {
      glm::dvec3 position = glm::dvec3(values[0], values[1], values[2]);
      glm::dvec3 color = glm::dvec3(values[3], values[4], values[5]);

      this->lights[this->numLights++] = new Light::PointLight(position, color, this->attenuation.top());
    }
  }
  else if (cmd == "pushTransform") {
    this->transforms.push(this->transforms.top());
  }
  else if (cmd == "popTransform") {
    if (this->transforms.size() <= 1) {
      std::cerr << "Stack has no elements.  Cannot Pop\n"; 
    } else {
      this->transforms.pop(); 
    }
  }
  else if (cmd == "translate") {
    if(readvals(lineStream, 3, values)) {
      Transform::rightmultiply(Transform::translate(values[0], values[1], values[2]),
                               this->transforms);
    }
  }
  else if (cmd == "rotate") {
    if(readvals(lineStream, 4, values)) {
      glm::dvec3 axis = glm::dvec3(values[0], values[1], values[2]);
      Transform::rightmultiply(Transform::rotate(values[3], axis), this->transforms);
    }
  }
  else if (cmd == "scale") {
    if(readvals(lineStream, 3, values)) {
      Transform::rightmultiply(Transform::scale(values[0], values[1], values[2]),
                               this->transforms);
    }
  }
}

Raytracer::Camera::Camera(glm::dvec3 eye, glm::dvec3 lookAt, glm::dvec3 up, double fovy) {
  this->eye = eye;
  this->lookAt = lookAt;
  this->up = up;
  this->fovy = fovy;
  this->computeOrthonormalBasis();
}

void Raytracer::Camera::computeOrthonormalBasis() {
  glm::dvec3 w = glm::normalize(this->eye - this->lookAt);
  glm::dvec3 u = glm::normalize(glm::cross(this->up, w));
  glm::dvec3 v = glm::cross(w, u);
  this->orthonormalBasis = glm::dmat3(u, v, w);
}

void Raytracer::Camera::print() {
  std::cout << "Camera eye:\n"; 
  std::cout << glm::to_string(this->eye);
  std::cout << "Camera lookAt:\n"; 
  std::cout << glm::to_string(this->lookAt);
  std::cout << "Camera up:\n"; 
  std::cout << glm::to_string(this->up);
  std::cout << "Camera fovy: " << this->fovy << "\n";
}
