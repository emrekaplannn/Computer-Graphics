#ifndef __COLOR_H__
#define __COLOR_H__

#include <iostream>

class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);
    //eklenti
    Color operator*(double scalar) const;
    Color operator/(double scalar);
    Color operator+(const Color &second);
    Color operator-(const Color &second);
    Color &operator=(const Color &second);
    Color operator+=(double second);


    friend std::ostream& operator<<(std::ostream& os, const Color& c);
};

#endif