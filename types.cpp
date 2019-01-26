#include "types.h"

Type::Cube::Cube(glm::dvec3 corners[Type::Cube::numCorners]) {
  for(int i = 0; i < Type::Cube::numCorners; i++) {
    this->corners[i] = corners[i];
  }
}

Type::AABB::AABB(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, Shape* shape) {
  this->xmin = xmin;
  this->xmax = xmax;
  this->ymin = ymin;
  this->ymax = ymax;
  this->zmin = zmin;
  this->zmax = zmax;
  this->shape = shape;
}

void Type::AABB::print() {
  std::cout << glm::to_string(glm::dvec3(this->xmin, this->ymin, this->zmin)) << "\n";
  std::cout << glm::to_string(glm::dvec3(this->xmax, this->ymax, this->zmax)) << "\n";
}

Type::Cube Type::AABB::toCube() {
  glm::dvec3 corners[Type::Cube::numCorners];

  corners[0] = glm::dvec3(this->xmin, this->ymin, this->zmin);
  corners[1] = glm::dvec3(this->xmin, this->ymin, this->zmax);
  corners[2] = glm::dvec3(this->xmin, this->ymax, this->zmin);
  corners[3] = glm::dvec3(this->xmin, this->ymax, this->zmax);
  corners[4] = glm::dvec3(this->xmax, this->ymin, this->zmin);
  corners[5] = glm::dvec3(this->xmax, this->ymin, this->zmax);
  corners[6] = glm::dvec3(this->xmax, this->ymax, this->zmin);
  corners[7] = glm::dvec3(this->xmax, this->ymax, this->zmax);

  return Type::Cube(corners);
}

bool Type::AABB::does_intersect(Type::Ray ray, Type::IntersectionInfo &hit) {
  double inverseDirectionX = 1.0 / ray.direction.x;
  double inverseDirectionY = 1.0 / ray.direction.y;
  double inverseDirectionZ = 1.0 / ray.direction.z;

	double t1 = (this->xmin - ray.origin.x) * inverseDirectionX;
  double t2 = (this->xmax - ray.origin.x) * inverseDirectionX;
  double t3 = (this->ymin - ray.origin.y) * inverseDirectionY;
  double t4 = (this->ymax - ray.origin.y) * inverseDirectionY;
  double t5 = (this->zmin - ray.origin.z) * inverseDirectionZ;
  double t6 = (this->zmax - ray.origin.z) * inverseDirectionZ;

  double tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
  double tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

  if (tmax < 0 || tmin > tmax) {
    return false;
  }

  hit.distance = tmin;
  return true;
}

Type::AABB* Type::Triangle::getAABB() {
  int numVertices = 3;
  glm::dvec3 transformedVertices[numVertices];

  transformedVertices[0] = Geometry::transformPoint(*this->v1, this->transform);
  transformedVertices[1] = Geometry::transformPoint(*this->v2, this->transform);
  transformedVertices[2] = Geometry::transformPoint(*this->v3, this->transform);

  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  double zmax = -std::numeric_limits<double>::infinity();
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double zmin = std::numeric_limits<double>::infinity();

  for(int i = 0; i < numVertices; i++) {
    xmin = std::min(transformedVertices[i].x, xmin);
    xmax = std::max(transformedVertices[i].x, xmax);
    ymin = std::min(transformedVertices[i].y, ymin);
    ymax = std::max(transformedVertices[i].y, ymax);
    zmin = std::min(transformedVertices[i].z, zmin);
    zmax = std::max(transformedVertices[i].z, zmax);
  }

  return new Type::AABB(xmin, ymin, zmin, xmax, ymax, zmax, this);
}

glm::dvec3 Type::Triangle::computeNormals() {
  glm::dvec3 planeVector1 = *this->v2 - *this->v1;
  glm::dvec3 planeVector2 = *this->v3 - *this->v1;
  this->normal = glm::normalize(glm::cross(planeVector1, planeVector2));
  this->transformedNormal = glm::normalize(
      Geometry::transformDirection(
        this->normal,
        this->transposedInverseTransform));
}

bool Type::Triangle::barycentric_test(Ray ray, double distance) {
    glm::dvec3 ray_segment = ray.getPoint(distance);
    glm::dvec3 v2 = ray_segment - *this->v1;
    glm::dvec3 v1 = *this->v2 - *this->v1;
    glm::dvec3 v0 = *this->v3 - *this->v1;
    double dot00 = glm::dot(v0, v0);
    double dot01 = glm::dot(v0, v1);
    double dot02 = glm::dot(v0, v2);
    double dot11 = glm::dot(v1, v1);
    double dot12 = glm::dot(v1, v2);
    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double beta = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double gamma = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return (beta >= 0) && (gamma >= 0) && (beta + gamma <= 1);
}

Type::Triangle::Triangle(glm::dvec3* vertex1, glm::dvec3* vertex2, glm::dvec3* vertex3, glm::dvec3 ambient, glm::dvec3 diffuse, glm::dvec3 emission, glm::dvec3 specular, double shininess, glm::dmat4 transform) {
  this->v1 = vertex1;
  this->v2 = vertex2;
  this->v3 = vertex3;
  this->material = Type::Material{ambient, diffuse, emission, specular, shininess};
  this->transform = transform;
  this->inverseTransform = glm::inverse(this->transform);
  this->transposedInverseTransform = glm::transpose(this->inverseTransform);
  this->computeNormals();
}

void Type::Triangle::print() {
  std::cout << "Triangle\n";
  std::cout << glm::to_string(*this->v1) << "\n";
  std::cout << glm::to_string(*this->v2) << "\n";
  std::cout << glm::to_string(*this->v3) << "\n";
  //printMaterial(this->material);
  //print(this->aabb);
}

bool Type::Triangle::does_intersect(Type::Ray ray, Type::IntersectionInfo &hit) {
  Type::Ray transformedRay = ray.transform(this->inverseTransform);

  double denominator = glm::dot(transformedRay.direction, this->normal);

  if(denominator == 0) {
    return false;
  }

  double vertexDotNormal = glm::dot(*this->v1, this->normal);
  double originDotNormal = glm::dot(transformedRay.origin, this->normal);
  double distance = (vertexDotNormal - originDotNormal) / denominator;

  glm::dvec3 point = Geometry::transformPoint(
      transformedRay.getPoint(distance),
      this->transform);

  if(distance < EPSILON || !this->barycentric_test(transformedRay, distance)) {
    return false;
  };

  hit = Type::IntersectionInfo(distance, point, this->transformedNormal, this);
  return true;
}

Type::Ray::Ray() {}

Type::Ray::Ray(glm::dvec3 o, glm::dvec3 d) {
  origin = o;
  direction = d;
}

Type::Ray Type::Ray::transform(glm::dmat4 transform) {
  glm::dvec4 homogenizedOrigin = glm::dvec4(this->origin, 1);
  glm::dvec4 homogenizedDirection = glm::dvec4(this->direction, 0);

  glm::dvec3 transformedOrigin = Geometry::dehomogenize(transform * homogenizedOrigin);
  glm::dvec3 transformedDirection = glm::dvec3(transform * homogenizedDirection);

  return Type::Ray(transformedOrigin, transformedDirection);
}

glm::dvec3 Type::Ray::getPoint(double distance) {
  return this->origin + this->direction * distance;
}

void Type::Ray::print()
{
  std::cout << "Ray\nOrigin\n";
  std::cout << glm::to_string(this->origin) << "\n";
  std::cout << "Direction\n";
  std::cout << glm::to_string(this->direction) << "\n";
}

Type::AABB* Type::Sphere::getAABB() {
  double xmin = this->center.x - this->radius;
  double xmax = this->center.x + this->radius;
  double ymin = this->center.y - this->radius;
  double ymax = this->center.y + this->radius;
  double zmin = this->center.z - this->radius;
  double zmax = this->center.z + this->radius;

  Type::AABB originalAABB(xmin, ymin, zmin, xmax, ymax, zmax, this);
  return originalAABB.toCube().transform(this->transform).getAABB(this);
}

Type::Sphere::Sphere(double x, double y, double z, double r, glm::dvec3 ambient, glm::dvec3 diffuse, glm::dvec3 emission, glm::dvec3 specular, double shininess, glm::dmat4 transform) {
  center = glm::dvec3(x, y, z);
  radius = r;
  material = Type::Material{ambient, diffuse, emission, specular, shininess};
  this->transform = transform;
  this->inverseTransform = glm::inverse(this->transform);
  this->transposedInverseTransform = glm::transpose(this->inverseTransform);
}

bool Type::Sphere::does_intersect(Type::Ray ray, Type::IntersectionInfo &hit) {
  Type::Ray transformedRay = ray.transform(this->inverseTransform);

  glm::dvec3 centerToRayOrigin = transformedRay.origin - this->center;
  double a = glm::dot(transformedRay.direction, transformedRay.direction);
  double b = 2.0f * glm::dot(transformedRay.direction, centerToRayOrigin);
  double c = glm::dot(centerToRayOrigin, centerToRayOrigin) - (this->radius * this->radius);

  double determinant = (b * b) - (4 * a * c);

  if (determinant < 0) {
    return false;
  }

  double determinant_square_root = sqrt(determinant);
  double root1 = (-b + determinant_square_root) / (2 * a);
  double root2 = (-b - determinant_square_root) / (2 * a);

  double distance;

  if (root1 > 0 && root2 > 0) {
    distance = root1 < root2 ? root1 : root2;
  } else if (root1 * root2 < 0) {
    distance = root1 > root2 ? root1 : root2;
  } else {
    return false;
  }

  glm::dvec3 intersectionPoint = transformedRay.getPoint(distance);

  glm::dvec3 point = Geometry::transformPoint(intersectionPoint, this->transform);

  glm::dvec3 centerToIntersection = glm::normalize(intersectionPoint - this->center);

  glm::dvec3 normal = glm::normalize(
      Geometry::transformDirection(
        centerToIntersection,
        this->transposedInverseTransform));

  hit = Type::IntersectionInfo(distance, point, normal, this);
  return true;
}

void Type::Sphere::print() {
  std::cout << "Sphere\n";
  std::cout << "Center:\n";
  std::cout << glm::to_string(this->center) << "\n";
  std::cout << "Radius: " << this->radius << "\n";
  //printMaterial(this->material);
  //print(this->aabb);
}


Type::AABB* Type::Cube::getAABB(Sphere* sphere) {
  double xmin, xmax, ymin, ymax, zmin, zmax;
  xmin = xmax = this->corners[0].x;
  ymin = ymax = this->corners[0].y;
  zmin = zmax = this->corners[0].z;

  for(int i = 1; i < Type::Cube::numCorners; i++) {
    xmin = std::min(this->corners[i].x, xmin);
    xmax = std::max(this->corners[i].x, xmax);
    ymin = std::min(this->corners[i].y, ymin);
    ymax = std::max(this->corners[i].y, ymax);
    zmin = std::min(this->corners[i].z, zmin);
    zmax = std::max(this->corners[i].z, zmax);
  }

  return new Type::AABB(xmin, ymin, zmin, xmax, ymax, zmax, (Shape*) sphere);
}

Type::Cube Type::Cube::transform(glm::dmat4 transforms) {
  glm::dvec3 corners[Type::Cube::numCorners];

  for(int i = 0; i < Type::Cube::numCorners; i++) {
    corners[i] = Geometry::transformPoint(this->corners[i], transforms);
  }

  return Type::Cube(corners);
}

Type::LineSegment** Type::Cube::getEdges()
{
  Type::LineSegment** edges = new LineSegment* [12];
  int edgeIndices[12][2] = {
    {0, 1},
    {0, 2},
    {0, 4},
    {3, 1},
    {3, 2},
    {3, 7},
    {6, 2},
    {6, 4},
    {6, 7},
    {5, 1},
    {5, 4},
    {5, 7}
  };

  for(int i = 0; i < 12; i++) {
    int start = edgeIndices[i][0];
    int end = edgeIndices[i][1]; 
    glm::dvec3 origin = this->corners[start];
    glm::dvec3 direction = this->corners[end] - this->corners[start];
    edges[i] = new LineSegment(origin, direction);
  }

  return edges;
}

Type::IntersectionInfo::IntersectionInfo(double distance, glm::dvec3 point, glm::dvec3 normal, Shape* shape) {
  this->distance = distance;
  this->point = point;
  this->normal = normal;
  this->shape = shape;
  this->displacedPoint = point + normal * EPSILON;
}

Type::IntersectionInfo::IntersectionInfo() {
  glm::dvec3 zero(0);
  this->distance = 0;
  this->point = zero;
  this->normal = zero;
  this->displacedPoint = zero;
  this->shape = NULL;
}

glm::dvec3 phongHighlight(
    glm::dvec3 &eye,
    Type::IntersectionInfo &hit,
    glm::dvec3 &directionToLight)
{
  glm::dvec3 directionToEye = glm::normalize(eye - hit.point);
  glm::dvec3 halfVec = glm::normalize(directionToLight + directionToEye);
  Type::Material material = hit.shape->material;

  double hDotN = glm::dot(halfVec, hit.normal);

  return material.specular * std::pow(std::max(hDotN, 0.0), material.shininess);
}

glm::dvec3 lambertianDiffusion(
    Type::IntersectionInfo &hit,
    glm::dvec3 &directionToLight)
{
  double nDotL = glm::dot(directionToLight, hit.normal);
  Type::Material material = hit.shape->material;

  return material.diffuse * std::max(nDotL, 0.0);
}

glm::dvec3 Type::Shape::computeColorAtPoint(
    glm::dvec3 &eye,
    Type::IntersectionInfo &hit,
    glm::dvec3 &directionToLight,
    glm::dvec3 &lightColor)
{
  glm::dvec3 specularComponent = phongHighlight(eye, hit, directionToLight);

  glm::dvec3 diffuseComponent = lambertianDiffusion(hit, directionToLight);

  return lightColor * (diffuseComponent + specularComponent);
}

Type::LineSegment::LineSegment() {}

Type::LineSegment::LineSegment(glm::dvec3 origin, glm::dvec3 direction) {
  this->ray = Type::Ray(origin, direction);
  glm::dvec3 white = glm::dvec3(1);
  this->material = Type::Material{white, white, white, white, 0.0};
}

void Type::LineSegment::print() {
  std::cout << "LineSegment\n";
  this->ray.print();
}

glm::dvec3 Type::LineSegment::computeColorAtPoint(
  glm::dvec3 &eye,
  Type::IntersectionInfo &hit,
  glm::dvec3 &directionToLight,
  glm::dvec3 &lightColor) 
{
  return glm::dvec3(1, 1, 1);
}

bool Type::LineSegment::does_intersect(Type::Ray ray, Type::IntersectionInfo &hit)
{
  glm::dvec3 o1 = this->ray.origin;
  glm::dvec3 d1 = this->ray.direction;

  glm::dvec3 o2 = ray.origin;
  glm::dvec3 d2 = ray.direction;

  glm::dvec3 originDiff = o2 - o1;
  glm::dvec3 directionCross = glm::cross(d1, d2);

  if(glm::all(glm::equal(directionCross, glm::dvec3(0)))) {
    return false;
  }

  double length = glm::length(directionCross);
  double denominator = glm::dot(length, length);

  double numerator = glm::dot(glm::cross(originDiff, d2), directionCross);
  double t1 = numerator / denominator;

  numerator = glm::dot(glm::cross(originDiff, d1), directionCross);
  double t2 = numerator / denominator;

  glm::dvec3 pointOnRay = ray.getPoint(t2);
  glm::dvec3 pointOnLineSegment = this->ray.getPoint(t1);
  double distanceBetweenPoints = glm::length(pointOnRay - pointOnLineSegment); 

  if(t2 >= 0 && t1 >= 0 && t1 <= 1.0 && distanceBetweenPoints < 0.009) {
    hit = Type::IntersectionInfo(t2, pointOnRay, glm::dvec3(0), this);
    return true;
  }

  return false;
}
