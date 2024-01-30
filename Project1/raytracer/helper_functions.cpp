
#include <iostream>
#include <cmath>
#include "helper_functions.h"

// firstly create vector operations systems

//default constructor
Vector_3d::Vector_3d() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
}

// default constructor (assign)
Vector_3d::Vector_3d(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}


//copy constructor
Vector_3d::Vector_3d(const Vector_3d &obj) {

    *this = obj;
}

// scalar multipication --- vektorun uzunlugunu degistiriyoruz.
Vector_3d
Vector_3d::operator*(float multiplier) {


    Vector_3d result;
    (result.x) = (this->x) * multiplier;
    (result.y) = (this->y) * multiplier;
    (result.z) = (this->z) * multiplier;
    return result;
}

//cross multipication ---- iki vektore dik ucuncu vektoru cikartiyoruz.
Vector_3d
Vector_3d::cross(const Vector_3d &new_vector) {
    Vector_3d result;
    result.x = (this->y) * new_vector.z - (this->z) * new_vector.y;
    result.y = (this->z) * new_vector.x - (this->x) * new_vector.z;
    result.z = (this->x) * new_vector.y - (this->y) * new_vector.x;
    return result;
}


// vector normalization --- we want to change the vector elements to <1
Vector_3d
Vector_3d::normalize() {
    Vector_3d result;
    float new_normalized_vector = (*this) * (*this);
    new_normalized_vector = 1.0 / sqrt(new_normalized_vector);
    result = (*this) * new_normalized_vector;
    return result;


}

//dot product of two vectors
float Vector_3d::operator*(const Vector_3d &new_vector) {

    return this->x * new_vector.x + this->y * new_vector.y + this->z * new_vector.z;

}

// addition of two vectors
Vector_3d
Vector_3d::operator+(const Vector_3d &new_vector) {

    Vector_3d result;
    (result.x) = (this->x) + new_vector.x;
    (result.y) = (this->y) + new_vector.y;
    (result.z) = (this->z) + new_vector.z;
    return result;
}

// deep copy constructor
Vector_3d &
Vector_3d::operator=(const Vector_3d &new_vector) {


    this->x = new_vector.x;
    this->y = new_vector.y;
    this->z = new_vector.z;
    return *this;

}