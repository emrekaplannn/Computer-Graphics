#ifndef _helperfunctions_h__
#define _helperfunctions_h__

#include <vector>
#include <cmath>
#include <iostream>

class Vector_3d;
class Ray;
class Sphere;
class Light;
class Triangle;
class Mesh;




class Vector_3d {
public:
    //vector cordinates

    float x, y, z;

    // vector constructors.

    Vector_3d();

    Vector_3d(float x, float y, float z);

    Vector_3d(const Vector_3d &obj);

    // vector calculations

    Vector_3d normalize();

    float len(const Vector_3d &point);

    Vector_3d cross(const Vector_3d &new_vector);
    
    // vector operators

    Vector_3d &operator=(const Vector_3d &new_vector);

    Vector_3d operator*(float multiplier);

    float operator*(const Vector_3d &new_vector);

    Vector_3d operator+(const Vector_3d &new_vector);


};

#endif

