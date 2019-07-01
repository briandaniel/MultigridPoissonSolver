/*%% Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Solves the laplace equation using finite differences on several grids   %
%                                                                         %
% Copyright (C) Brian D Hong, 2012                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef HEADER_FILE
#define HEADER_FILE

// Includes
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex>
#include <iomanip> 
#include <fstream>
using namespace std;

// Complex variable type definition "CX"
typedef std::complex<double> CX;

// Include class definitions
#include "finiteDifference.hpp" 
#include "mathFunctions.hpp"

// Use open MP for parallel computation?
// #define USE_OPENMP
#ifdef USE_OPENMP
	#include <omp.h>
	#define OMP_NUM_THREADS 16
#endif


// Standard constant definitions
#define PI	   3.14159265
#define EPSNOT 8.854187817620e-12
#define MUNOT  1.25663706e-6
#define CNOT   2.99792458e8
#define ETANOT 3.767303132465403e02
#define PLANCK 6.62606957e-34
#define E      1.60217646e-19
#define ME     9.10938188e-31
#define MP     1.66053892e-27
#define HBAR   1.0545717253e-34  
#define JEV    6.24150974e18


#endif // HEADER_FILE
