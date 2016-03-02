#include "RealTimeMath.h"
#include <math.h>
#include <iostream>
using namespace std;
RealTimeMath::RealTimeMath(){
    powFast = new PowFast(18);
    }

RealTimeMath::~RealTimeMath(){
    delete powFast;
    }

double RealTimeMath::fastEXP(double x){
	if(x > 88.5 || x < -87){ // Overflow prevention due to powFast library's use of float
		if( x < -746 ) // Close to 0 shortcut: If x < -746, math.h exp(x) function returns 0 due to truncation
			return 0;

		if( x > 710 ){ // Infinity shortcut: If x > 710, math.h exp(x) returns infinity
			return INFINITY;
		}

		// Standard overflow prevention
		// Due to powFast limitations from using floats,
		// if x is too large or small, problem is broken into smaller pieces
		// e.g. fastEXP(100) ==> e^180 = e^60 * e^60 * e^60
		powCount = 1;
		powAns = 1;
		powExp = x;

		while(powExp > 88.5 || powExp < -87){
			powCount++;
			powExp = x / powCount;
		}

		for(n = 1; n <= powCount; n++){
			powAns = powAns * powFast->e(powExp);
		}

		return powAns;
		
	} // end overflow prevention
		
	else // call powFast e^x function
		return powFast->e(x);
}

// POWER OPTIMIZATION
// x^y = e^(y*ln(x))

double RealTimeMath::fastPOW(double x, double y){
    return fastEXP(y * log(x)); // Pow function using fastEXP
}
