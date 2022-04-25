#include"stiff.cc"
static const real k = 2e+5;

int main(int argc, char** argv)
{
Stiff s;
std::vector<real> Y;
std::vector<std::string> str = s.ReadParams("Input.txt");

int flag = stoi(str[0]);
real x0 , xf, ic, dt;
x0 = stof(str[1]);
xf = stof(str[2]);
dt = stof(str[3]);
ic = stof(str[4]);

auto f = [&k](real y, real x, real dt){ return k*( -y + exp(-x) ) - exp(-x) ; } ; 
auto g = [&k](real y, real x, real dt){ return ( y + k*dt*exp(-x) - dt*exp(-x) )/(1 + k*dt) ;};


Stiff stf(x0, xf, dt, ic);

if(flag == 1)
{
	Y = stf.Euler(f,flag);
	stf.WriteSol(Y, "solution.dat");
}
else if(flag == 2)
{
	Y = stf.Euler(g,flag);
	stf.WriteSol(Y, "solution.dat");
}
else
{
	std::cout<<"Incorrect flag Value \n Enter 1 for Explicit Euler or 2 for Implicit Euler\n"<<std::endl;
}


return EXIT_SUCCESS;
}
