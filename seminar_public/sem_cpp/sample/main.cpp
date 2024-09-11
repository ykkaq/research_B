#include <iostream>
#include "interval.hpp"

template <typename _T>
_T add( _T a, _T b ){
    return a + b;
}

int main(void){
    interval< float > a, b, c;
    a.setting(1.0f, 2.0f);
    b.setting(2.0f, 3.0f);

    c =  b - a + b + b + a - b;
    // c = a.add(b);
    c.disp();
    c = sqrt(c);
    c.disp();
}
