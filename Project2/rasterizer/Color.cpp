#include "Color.h"
#include <iostream>
#include <iomanip>

using namespace std;

Color::Color() {}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}
// helper fonksiyonlar
Color
Color::operator*(double scalar) const{

    Color result;
    (result.r) = (this->r) * scalar;
    (result.g) = (this->g) * scalar;
    (result.b) = (this->b) * scalar;
    return result;
}

Color
Color::operator+=(double scalar){

    Color result;
    (result.r) = (this->r) + scalar;
    (result.g) = (this->g) * scalar;
    (result.b) = (this->b) * scalar;
    return result;
}

Color
Color::operator/(double scalar) {

    Color result;
    (result.r) = (this->r) / scalar;
    (result.g) = (this->g) / scalar;
    (result.b) = (this->b) / scalar;
    return result;
}

Color
Color::operator+(const Color &second) {

    Color result;
    (result.r) = (this->r) + second.r;
    (result.g) = (this->g) + second.g;
    (result.b) = (this->b) + second.b;
    return result;
}
Color
Color::operator-(const Color &second) {

    Color result;
    (result.r) = (this->r) - second.r;
    (result.g) = (this->g) - second.g;
    (result.b) = (this->b) - second.b;
    return result;
}
Color &
Color::operator=(const Color &second) {


    this->r = second.r;
    this->g = second.g;
    this->b = second.b;
    return *this;

}


ostream& operator<<(ostream& os, const Color& c)
{
    os << fixed << setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}
