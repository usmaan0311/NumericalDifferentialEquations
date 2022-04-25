#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<cassert>
#include<cmath>
#include<functional>

typedef double real;
class Stiff{

real x0;
real xf;
real h;
real ic;
long int n;
std::vector<real> X, Y;

	public:
		Stiff();
		Stiff(real, real, real, real);
		std::vector<std::string> ReadParams(std::string);
		template<typename T>
		std::vector<real> Euler(T,int); 
//		std::vector<real> Euler(std::function<real(real,real,real)>, std::function<real(real,real)>,int);
		void WriteSol(std::vector<real> , std::string);

};
