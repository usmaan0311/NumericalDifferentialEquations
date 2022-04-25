#include<iostream>
#include"sms.cc"


//////////////////////////////////////////////////////////////////////////////////////////////////
//												//
//												//
//				m*(d2x/dt2) + c*(dx/dt) + k*x = 0				//
//												//	
//					dx/dt = z	....(1)					//
//												//
//				dz/dt = -( c*z + k*x )/m	.....(2)			//
//												//
//												//		
/////////////////////////////////////////////////////////////////////////////////////////////////


static const real m=20.0, k=20.0;


int main(int argc, char* argv[])
{
	Sms s;
	std::vector<std::string>	ip = s.InpParams("Input.txt");

	real c = stof(ip[0]);
	real x0 = stof(ip[1]);
	real v0 = stof(ip[2]);
	real dt = stof(ip[3]);
	real ti=0, tf=15;
	Sms sol(x0, v0, ti, tf, dt);
	auto f = [](real t, real x, real z){return z;};
	auto g = [c](real t, real x, real z){return -(c*z + k*x)/m ;};

	std::vector<real> Xc = sol.RK4(f,g)[0];
	sol.WriteSol(Xc, "solution.dat");

return EXIT_SUCCESS;
}

