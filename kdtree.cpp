#include "kdtree.h"

Type::AABB* getEnclosingAABB(Type::AABB** aabbs, int start, int end) {
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  double zmax = -std::numeric_limits<double>::infinity();
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double zmin = std::numeric_limits<double>::infinity();

  for(int i = start; i < end; i++) {
    xmin = std::min(aabbs[i]->xmin, xmin);
    xmax = std::max(aabbs[i]->xmax, xmax);
    ymin = std::min(aabbs[i]->ymin, ymin);
    ymax = std::max(aabbs[i]->ymax, ymax);
    zmin = std::min(aabbs[i]->zmin, zmin);
    zmax = std::max(aabbs[i]->zmax, zmax);
  }

  return new Type::AABB(xmin, ymin, zmin, xmax, ymax, zmax, NULL);
}

void swap(Type::AABB** aabbs, int i, int j) {
  Type::AABB* temp = aabbs[i];
  aabbs[i] = aabbs[j];
  aabbs[j] = temp;
}

void sortByX(Type::AABB** aabbs, int start, int end) {
  for(int i = start; i < end; i++) {
    for(int j = i; j < end; j++) {
      if(aabbs[i]->xmin > aabbs[j]->xmin) {
        swap(aabbs, i, j);
      }
    }
  }
}

void sortByY(Type::AABB** aabbs, int start, int end) {
  for(int i = start; i < end; i++) {
    for(int j = i; j < end; j++) {
      if(aabbs[i]->ymin > aabbs[j]->ymin) {
        swap(aabbs, i, j);
      }
    }
  }
}

void sortByZ(Type::AABB** aabbs, int start, int end) {
  for(int i = start; i < end; i++) {
    for(int j = i; j < end; j++) {
      if(aabbs[i]->zmin > aabbs[j]->zmin) {
        swap(aabbs, i, j);
      }
    }
  }
}

Type::AABB* KDTree::buildKdTree(Type::AABB** aabbs, int start, int end, int axis) {
  if(end - start == 1) {
    return aabbs[start];
  }

  Type::AABB* parent = getEnclosingAABB(aabbs, start, end);

  switch(axis) {
    case 0:
      sortByX(aabbs, start, end);
      break;
    case 1:
      sortByY(aabbs, start, end);
      break;
    case 2:
      sortByZ(aabbs, start, end);
      break;
  }

  int nextAxis = (axis + 1) % 3;
  int half = (start + end) / 2;

  parent->children[0] = buildKdTree(aabbs, start, half, nextAxis);
  parent->children[1] = buildKdTree(aabbs, half, end, nextAxis);

  return parent;
}
