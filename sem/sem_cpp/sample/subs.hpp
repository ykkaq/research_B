#include <iostream>
#include <cfenv>
#include <cmath>
#include "add.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template< typename _T>
struct isubs : public iadd< _T >{
    void subs_aminusb( const isubs< _T >& b ){
        volatile _T ainf = this->inf, asup = this->sup;
        volatile _T binf = b.inf, bsup = b.sup;
        volatile _T rinf, rsup;

        fesetround(FE_DOWNWARD);
        rinf = ainf - bsup;
        
        fesetround(FE_UPWARD);
        rsup = asup - binf;

        fesetround(FE_TONEAREST);
        this->sup = rsup;
        this->inf = rinf;
    }

    void subs_bminusa( const isubs< _T >& b ){
        volatile _T ainf = this->inf, asup = this->sup;
        volatile _T binf = b.inf, bsup = b.sup;
        volatile _T rinf, rsup;

        fesetround(FE_DOWNWARD);
        rinf = binf - asup;
        
        fesetround(FE_UPWARD);
        rsup = bsup - ainf;

        fesetround(FE_TONEAREST);
        this->sup = rsup;
        this->inf = rinf;
    }
};