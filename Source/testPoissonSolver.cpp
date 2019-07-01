/*%% Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Solves the laplace equation using finite differences on several grids   %
%                                                                         %
% Copyright (C) Brian D Hong, 2012                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "testPoissonSolver.hpp"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int main () {


	int Nx, Ny, Nz, N3;
	double Lx, Ly, Lz;
	double dx, dy, dz;
	string boundary;

	boundary = "zero";

	double K = 7.0;
	
	Nx = (int) pow(2.0, K) - 1;
	Ny = (int) pow(2.0, K) - 1;
	Nz = (int) pow(2.0, K) - 1;

	Nx = 127;
	Ny = 127;
	Nz = 127;
	
	N3 = Nx*Ny*Nz;

	Lx = 6;
	Ly = 6;
	Lz = 6;

	dx = Lx/(Nx-1);
	dy = Ly/(Ny-1);
	dz = Lz/(Nz-1);

	double* U = new double [N3];
	double* Uexact = new double [ N3];
	double* f = new double [N3];
	double* dU = new double [N3];
	double* guess = new double [N3];

	for (int k = 0; k < Nz; k++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
				double z = dz*k;
				double y = dy*j;
				double x = dx*i;

		
				// A data set
				/*
				Uexact[ k + j*Nz + i*Ny*Nz ] = exp ( -pow(x -3, 2.0) -pow(y -3, 2.0) -pow(z - 3, 2.0) );
				f[ k + j*Nz + i*Ny*Nz ] = 2.0* exp ( -pow(x -3, 2.0) -pow(y -3, 2.0) -pow(z - 3, 2.0) ) * 
					(2.0 * pow(x,2.0) - 12.0*x + 2.0 * pow(y,2.0) - 12.0*y + 2.0 * pow(z,2.0) - 12.0*z +51.0);
				guess[ k + j*Nz + i*Ny*Nz ] = timeRand();
				*/

				// Another data set
					Uexact[ k + j*Nz + i*Ny*Nz ] = sin(6*x)*exp ( -pow(x -3, 2.0) -pow(y -3, 2.0) -pow(z - 3, 2.0) );

				f[ k + j*Nz + i*Ny*Nz ] =  exp ( -pow(x -3, 2.0) -pow(y -3, 2.0) -pow(z - 3, 2.0) ) * 
					( 2.0* sin(6*x) * (2.0 * pow(x,2.0) - 12.0*x + 2.0 * pow(y,2.0) - 12.0*y + 2.0 * pow(z,2.0) - 12.0*z +33.0) - 24*(x-3)*cos(6*x) );

				guess[ k + j*Nz + i*Ny*Nz ] = timeRand();



			}
		}
	}
	
	Derivative D2(2.0, 2.0 , boundary, Nx, Ny, Nz, Lx, Ly, Lz);
	D2.deriv(Uexact, dU);
	PoissonMG pSolve( boundary, Nx, Ny, Nz, Lx, Ly, Lz);

	double elapTime;
	clock_t beginT, endT;
	beginT = clock();


	pSolve.PoissonSolve(f, guess, U);

	endT = clock();

	elapTime = ((double)(endT - beginT))/CLOCKS_PER_SEC;
	cout << elapTime << "[s]" << endl;

	////////////////////////////////////////////////////////

	int i = 0;

	N3 = Nx*Ny*Nz;
	std::ofstream intData("intData.bin", ios::out | ios::binary);
	std::ofstream data("data.bin", ios::out | ios::binary);

	for (int k = 0; k < N3; k++)
		data.write((char *) & U[k], sizeof U[k]);

	for (int k = 0; k < N3; k++)
		data.write((char *) & f[k], sizeof f[k]);

	for (int k = 0; k < N3; k++)
		data.write((char *) & Uexact[k], sizeof Uexact[k]);

	for (int k = 0; k < N3; k++)
		data.write((char *) & guess[k], sizeof guess[k]);

	for (int k = 0; k < N3; k++)
		data.write((char *) & dU[k], sizeof dU[k]);


	intData.write((char *) & N3, sizeof N3);
	intData.write((char *) & Nx, sizeof Nx);
	intData.write((char *) & Ny, sizeof Ny);
	intData.write((char *) & Nz, sizeof Nz);

	intData.close();
	data.close();

	////////////////////////////////////////////////////////

	cout << "Program finished." << endl;

	return (0);
}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


