/*
 * mathfunctions.hpp
 *
 *  Created on: Jul 1, 2019
 *      Author: brian
 */

#ifndef MATHFUNCTIONS_HPP_
#define MATHFUNCTIONS_HPP_

#include <math.h>
#include <complex>

// Complex variable type definition "CX"
typedef std::complex<double> CX;

/* mathFunctions.cpp declarations
OUTPUT                      FUNCTION		( INPUT )																						*/
double						round			( double d )																					;
double						trapz			( double *F, int N, double dz )																	;
CX							CXtrapz			( CX *F, int N, double dz )																		;
double*						linspace		( double a, double b, int N )																	;
double						dot				( double* U, double* V, int N )																	;
double						norm2			( CX* U, int N )																				;
double						mod				(double x, double y )																			;
void						sqNormalize		( CX** PSI, int N, int Npsi, double dz )														;
void 						CXtridisolve	( CX* PSI, CX* a, CX* b, CX* c, CX* d, CX* e, CX* f, int N )									;
void						tridisolve		( double* PSI, double* a, double* b, double* c, double* d, double* e, double* f, int N )		;
void						fftConv			( double* fgConv, double* f, double* g, int N, int xo )											;
void						fftShift		( double* shiftedg, double* g, int N )															;
void						derivative		( double* dF, double* F, double dz, int N )														;
void						spline			( double* x, double* y, double* xx, double* yy, int N, int M )									;
void						CXderivative2	( CX* d2F, CX* F, double dz, int N )															;
void						derivative2		( double* d2F, double* F, double dz, int N )													;
void						xySort			( double** x, double* y,  int M, int N )														;
double						timeRand		()																								;




#endif /* MATHFUNCTIONS_HPP_ */
