import math
import numpy as np

def rotate(degrees, axis):
    angleInRadians = math.radians(degrees)
    cosine = math.cos(angleInRadians)
    sine = math.sin(angleInRadians)

    identity = np.identity(3)

    x = axis[0]
    y = axis[1]
    z = axis[2]

    dual = np.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]]);

    rotationMatrix = cosine * identity + (1 - cosine) * np.outer(axis, axis) + sine * dual;

    homogenizedRotation = np.r_[np.c_[rotationMatrix, [0, 0, 0]], [[0, 0, 0, 1]]];

    return homogenizedRotation


camera = np.array([6, 2, 4, 1])
camera_demo = np.array([10, 10, 14, 1])
keytransform = rotate(-45, np.array([0, 1, 0])) * rotate(-45, np.array([1, 0, 0]))
filltransform = rotate(30, np.array([0, 1, 0])) * rotate(-20, np.array([1, 0, 0]))
print("Key light position: ", camera_demo * keytransform)
print("Fill light position: ", camera_demo * filltransform)
