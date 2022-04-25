#pragma once
#include<fstream>
#include<vector>
#include<string>
#include<cassert>
#include<cmath>
#include<functional>

typedef double real;

class Sms
{
real x0;
real v0;
real ti;
real tf;
real h;
int n;
std::vector<real> T,Y;

	public:
		Sms();
		Sms(real, real, real, real, real);
		std::vector<std::string> InpParams(std::string);
		std::vector<std::vector<real>> RK4(std::function<real(real, real, real)>, std::function<real(real, real, real)> );
		void WriteSol(std::vector<real>, std::string);
		void WriteErr(std::vector<real>, std::string);

};
