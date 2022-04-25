#include"sms.h"

Sms::Sms()
{}

Sms::Sms(real x0, real v0, real ti, real tf, real h)
	:x0(x0), v0(v0), ti(ti), tf(tf), h(h)
{
n = static_cast<int>( (tf - ti)/h );
for(int i=0; i<n;i++)
{
Y.emplace_back(0);
T.emplace_back(i*h);
}

}

std::vector<std::string> Sms::InpParams(std::string fname)
{
	std::ifstream fin;
	fin.open(fname,std::ios::in);

	std::vector<std::string> SVec;
	std::string str;
	while(std::getline(fin,str))
	{
		std::getline(fin,str);
		if((str[0] !='/' && str[0] !='*')&&(str[1] !='/' && str[1] !='*'))
		{
			SVec.emplace_back(str);

		}



	}

return SVec;
}


std::vector<std::vector<real>> Sms::RK4(std::function<real(real t, real x, real z)> f, std::function<real(real t, real x, real z)> g )
{
std::vector<real> Z = Y;
std::vector<std::vector<real>> SolVec;
// Initial conditions

Y[0]=x0 ; // y(x=0)
Z[0]=v0 ; // y'(x=0)

// RK4 implementation 

for(int i=0; i<n-1; i++)
{
real k1 = h*f(T[i],Y[i],Z[i]); real l1 = h*g(T[i],Y[i],Z[i]);
real k2 = h*f(T[i] + h/2, Y[i] + k1/2, Z[i] + l1/2) ; real l2 = h*g(T[i] + h/2, Y[i] + k1/2, Z[i] + l1/2);
real k3 = h*f(T[i] + h/2, Y[i] + k2/2, Z[i] + l2/2) ; real l3 = h*g(T[i] + h/2, Y[i] + k2/2, Z[i] + l2/2);
real k4 = h*f(T[i] + h, Y[i] + k3, Z[i] + l3) ; real l4 = h*g(T[i] + h, Y[i] + k3, Z[i] + l3);

Y[i+1] = Y[i] + (k1 + 2*(k2 + k3) + k4)/6 ; Z[i+1] = Z[i] + (l1 + 2*(l2 + l3) + l4)/6 ;


}
SolVec.emplace_back(Y);
SolVec.emplace_back(Z);

return SolVec;

}

void Sms::WriteSol(std::vector<real> Y, std::string Yname)
{
	std::ofstream YWrite(Yname);
	YWrite.precision(10);
	YWrite.setf(std::ios::scientific);
	assert(YWrite.is_open());
	for(int i=0;i<n;i++)
	{
		YWrite<<T[i]<<" "<<Y[i]<<"\n";
	
	}
	YWrite.close();


}

void Sms::WriteErr(std::vector<real> E, std::string Ename)
{
	std::ofstream EWrite(Ename);
	EWrite.precision(10);
	EWrite.setf(std::ios::scientific);
	assert(EWrite.is_open());
	for(int i=0;i<n;i++)
	{
		EWrite<<T[i]<<" "<<Y[i]<<"\n";
	
	}
	EWrite.close();


}

