#include<iostream>
#include<cmath>
#include"sol2.cc"

static real k=1.1;

real ya(real t);
real f(real t, real y);

int main(int argc, char* argv[])
{
std::vector<real> Yc,Err;

// default constructor for reading input from file
Sol2 sol;

//reading parameters from the input file
std::vector<std::string> params = sol.ReadParams("Input.txt");


int flag = stoi(params[0]);
real x0 = stof(params[1]);
real xf = stof(params[2]);
real ic = stof(params[3]);
real dt = stof(params[4]);

Sol2 sol1(ic,dt,x0,xf);


// lambda function, analytical solution for the differential equation
/*auto ya=[k](real t){return exp(pow(t,3)/3 - k*t );};
auto f=[k](real t, real y){return y*pow(t,2) - k*y;};
*/


if(flag==1)
{
Yc = sol1.Euler(f);
Err = sol1.AvgErr(ya,f,flag)[1];
sol1.WriteSol(Yc,"solution.dat",flag);
sol1.WriteErr(Err,"Error.dat",flag);
}
else if(flag==2)
{
Yc = sol1.AB2(f);
Err = sol1.AvgErr(ya,f,flag)[1];
sol1.WriteSol(Yc,"solution.dat",flag);
sol1.WriteErr(Err,"Error.dat",flag);
}
else
{
Yc = sol1.RK4(f);
Err = sol1.AvgErr(ya,f,flag)[1];
sol1.WriteSol(Yc,"solution.dat",flag);
sol1.WriteErr(Err,"Error.dat",flag);
}


return EXIT_SUCCESS;
}



real ya(real t)
{
	return exp( pow(t,3)/3 - k*t );
}
real f(real t, real y)
{
return y*pow(t,2) - k*y;
}
