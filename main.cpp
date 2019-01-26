#include "raytracer.h"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: raytrace scenefile\n"; 
    exit(1); 
  }

  Raytracer::Scene scene(argv[1]);

  //scene.print();

  scene.draw();

  return 0;
}
