#pragma once

#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include<cassert>
#include<fstream>
#include<thread>
#include<algorithm>

typedef double real;

class Laplace{

int Nx,Ny;
int NST, TNX, TNY, PNX1, PNX2, PNY1, PNY2;
real dx,dy;
real xspan, yspan;
real DBC1, DBC2;
std::vector<real> x,y;
std::vector<std::vector<real>> T;
	public:
	Laplace();
	Laplace(int,int,int,real,real, real,real);
	real l2_norm(std::vector<std::vector<real>>, std::vector<std::vector<real>>);
	std::vector<std::vector<real>> Jacobi(real, int, std::vector<int>*, std::vector<real>*);
	std::vector<std::vector<real>> GaussSeidel(real, int, std::vector<int>*, std::vector<real>*);
	std::vector<std::vector<real>> SOR(real, int, std::vector<int>*, std::vector<real>*, real);
	void PrintSol(std::vector<std::vector<real>>, std::vector<real>, std::vector<int>, std::string);
	void WriteSol(std::vector<std::vector<real>>, std::string);
	void WriteXY(std::string, std::string);
	void WriteResidual(std::vector<real>, std::string);
	void WriteIter(std::vector<int>, std::string);

};

