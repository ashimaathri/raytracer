CC = g++
CFLAGS = -g -fopenmp
LDFLAGS = -lfreeimage
INCFLAGS = -I./glm-0.9.7.1

RM = rm -f
all: raytrace
raytrace: main.o types.o geometry.o transforms.o light.o raytracer.o kdtree.o
	$(CC) $(CFLAGS) -o raytrace main.o types.o geometry.o transforms.o light.o raytracer.o kdtree.o $(LDFLAGS)
main.o: main.cpp raytracer.h
	$(CC) $(CFLAGS) -c main.cpp
types.o: types.cpp types.h geometry.h light.h transforms.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c types.cpp
geometry.o: geometry.cpp geometry.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c geometry.cpp
transforms.o: transforms.cpp transforms.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c transforms.cpp
light.o: light.cpp light.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c light.cpp
raytracer.o: raytracer.cpp raytracer.h types.h kdtree.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c raytracer.cpp
kdtree.o: kdtree.cpp kdtree.h types.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c kdtree.cpp
clean: 
	$(RM) *.o raytrace *.png
