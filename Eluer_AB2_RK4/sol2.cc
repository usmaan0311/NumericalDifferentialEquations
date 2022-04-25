#include<string>
#include<vector>
#include<fstream>
#include<cassert>
#include"sol2.h"

Sol2::Sol2()
{}

Sol2::Sol2(real ic, real h,real x0, real xf)
	:ic(ic), h(h), x0(x0), xf(xf)
{

n = static_cast<int>( (xf - x0)/h) ;

for(int i=0;i<n;i++)
{
Y.emplace_back(0);
T.emplace_back(i*h);
}

}

std::vector<std::string> Sol2::ReadParams(std::string fname)
{
	std::ifstream fin;
	fin.open(fname,std::ios::in);

	std::vector<std::string> Rvec;
	std::string strg;

	while(std::getline(fin,strg))
	{
		std::getline(fin,strg);
	if( (strg[0] !='/' && strg[0] != '*') && (strg[1] != '/' && strg[1] != '*') )
		Rvec.push_back(strg);
	
	
	}

return Rvec;

}





std::vector<real> Sol2::Euler(real (*f)(real t, real y))
{
	std::vector<real> Ye = Y;
	Ye[0] = ic;
	for(int i=0; i<n;i++)
	{
	Ye[i+1]=Ye[i] + h*( f(T[i],Ye[i]) );
	}

	return Ye;

}

std::vector<real> Sol2::AB2(real (*f)(real t, real y))
{

	std::vector<real> Yeu,Yab=Y;
	Yeu=Euler(f);

	Yab[0] = ic;
	Yab[1] = Yeu[1];
	for(int i=1; i<n;i++)
	{
	Yab[i+1]=Yab[i] + h*( 1.5*f(T[i],Yab[i]) - 0.5*f(T[i-1],Yab[i-1]) );
	
	}
	return Yab;

}

std::vector<real> Sol2::RK4(real (*f)(real t,real y))
{
	real k1, k2, k3, k4;
	std::vector<real> Yrk=Y;
	Yrk[0] = ic;
	for(int i=0; i<n;i++)
	{
	k1 = h*f(T[i],Yrk[i]);
	k2 = h*f(T[i] + h/2, Yrk[i] + k1/2);
	k3 = h*f(T[i] + h/2, Yrk[i] + k2/2);
	k4 = h*f(T[i] + h, Yrk[i] + k3);
	Yrk[i+1] = Yrk[i] + (k1 + 2*(k2 + k3) + k4)/6;
	
	}
	return Yrk;

}


std::vector<std::vector<real>> Sol2::AvgErr(real (*fa)(real t), real (*f)(real t, real y), int flag)
{

	std::vector<std::vector<real>> ErrVec;
	std::vector<real> MaxEr, AvgEr,Yc;
	MaxEr=Y;
	AvgEr=Y;

	if(flag ==1)
	Yc = Euler(f);
	else if(flag==2)
	Yc = AB2(f);
	else
	Yc = RK4(f);
	
	for(int i=0;i<n;i++)
	{
	MaxEr[i] = fabs( fa(T[i]) - Yc[i] );
	AvgEr[i] = MaxEr[i]/n;
	}

	ErrVec.emplace_back(MaxEr);
	ErrVec.emplace_back(AvgEr);
return ErrVec;

}




void Sol2::WriteSol(std::vector<real>& Y, std::string Yname, int flag)
{
// Writing Solution
	std::ofstream WriteY(Yname);
	WriteY.precision(10);
	WriteY.setf(std::ios::scientific);
	assert(WriteY.is_open());
	if(flag==1)
		WriteY<<"# Explicit Euler"<<"\n";
	else if(flag==2)
		WriteY<<"# Second Order Adams-Bashforth"<<"\n";
	else
		WriteY<<"# Fourth-order Runge-Kutta"<<"\n";
	for(int i=0; i<n;i++)
	{
		WriteY<<T[i]<<" "<<Y[i]<<"\n";
	
	}
	WriteY.close();
}


void Sol2::WriteErr(std::vector<real>& Err, std::string Ename, int flag)
{
// Writing Errors


	std::ofstream WriteE(Ename);
	WriteE.precision(10);
	WriteE.setf(std::ios::scientific);
	assert(WriteE.is_open());
	if(flag==1)
		WriteE<<"# Explicit Euler"<<"\n";
	else if(flag==2)
		WriteE<<"# Second Order Adams-Bashforth"<<"\n";
	else
		WriteE<<"# Fourth-order Runge-Kutta"<<"\n";
	for(int i=0; i<n;i++)
	{
		WriteE<<T[i]<<" "<<Err[i]<<"\n";
	
	}
	WriteE.close();



}
