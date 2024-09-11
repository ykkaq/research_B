#include <iostream>
#include <cfenv>
#include <cmath>
#include "subs.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template <typename _T>
class interval : protected isubs< _T >{
public:
    void setting( const _T& i, const _T& s ){
        this->inf = i;
        this->sup = s;
    }

    void disp(void)const{
         std::cout << "[" << this->inf << "," << this->sup << "]" << std::endl; 
    }

    friend interval< _T > operator+( const interval< _T >& a, const interval< _T >& b ){
        interval< _T > res = a;
        res.add(b);
        return res;
    }

    friend interval< _T > operator+( interval< _T >&& a, const interval< _T >& b ){
        a.add(b);
        return std::move(a);
    }

    friend interval< _T > operator+( const interval< _T >& a, interval< _T >&& b ){
        b.add(a);
        return std::move(b);
    }

    friend interval< _T > operator+( interval< _T >&& a, interval< _T >&& b ){
        a.add(b);
        return std::move(a);
    }

    friend interval< _T > operator-( const interval< _T >& a, const interval< _T >& b ){
        interval< _T > res = a;
        res.subs_aminusb(b);
        return res;
    }

    friend interval< _T > operator-( interval< _T >&& a, const interval< _T >& b ){
        a.subs_aminusb(b);
        return std::move(a);
    }

    friend interval< _T > operator-( const interval< _T >& a, interval< _T >&& b ){
        b.subs_bminusa(a);
        return std::move(b);
    }

    friend interval< _T > operator-( interval< _T >&& a, interval< _T >&& b ){
        a.subs_aminusb(b);
        return std::move(a);
    }

    template <typename _TM> 
    friend std::enable_if< std::is_constructible< _T, _TM >::value, interval< _T > > operator+(const _TM& a, const interval< _T >& b ){
        _T aa = _T(a);
        aa.add(b);
        return aa;
    }

    /*friend interval< _T > operator-( const interval< _T >& a, const interval< _T >& b ){
        interval< _T > res = a;
        res.subs(b);
        return res;
    }*/

// 野良関数をclass内に書くにはfriend
//     → メソッド関数ではない!! * a.xxxx(  )
    friend interval< _T > sqrt( const interval< _T >& a ){
        interval< _T > res;
        volatile _T asup = a.sup, ainf = a.inf;
        volatile _T rsup, rinf;
        fesetround(FE_UPWARD);
        rsup = std::sqrt(asup);
        fesetround(FE_DOWNWARD);
        rinf = std::sqrt(ainf);
        fesetround(FE_TONEAREST);
        res.sup = rsup;
        res.inf = rinf;
        return res;
    }
};


struct F {
    friend double operator () ( double a ){
        return a*a;
    }
};

F f;
f(3)