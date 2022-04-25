#include"stiff.h"

Stiff::Stiff()
{}

Stiff::Stiff(real x0, real xf, real h, real ic)
	: x0(x0), xf(xf), h(h), ic(ic)
{
n = (int) (xf - x0)/h ;

for(int i=0; i<n; i++)
{
Y.emplace_back(0) ;
X.emplace_back(i*h) ;
}

}

std::vector<std::string> Stiff::ReadParams(std::string fname)
{
	std::ifstream fin;
	fin.open(fname, std::ios::in);
	std::vector<std::string> SVec;
	std::string str;

	while(std::getline(fin,str))
	{
		std::getline(fin,str);
		if( (str[0] != '/' && str[0] !='*')&& (str[1] != '/' && str[1] != '*') )
			SVec.emplace_back(str);

	}
	
return SVec;

}

template<typename T> std::vector<real> Stiff::Euler(T f, int flag)
{
std::cout<<"flag\t:"<<flag<<std::endl;
Y[0] = 1;
if(flag==1)
{
for(int i=0;i<n-1;i++)
{
Y[i+1] = Y[i] +  h*f(Y[i],X[i],h) ;
}
}
else
{
	for(int i=0;i<n-1;i++)
	{	Y[i+1] = f(Y[i],X[i],h); 
	}

}

 return Y;

}


void Stiff::WriteSol(std::vector<real> Y, std::string fname)
{
	std::ofstream fout(fname);
	fout.precision(10);
	fout.setf(std::ios::scientific);
	assert(fout.is_open());
	for(int i=0; i<n; i++)
		{
		fout<<X[i]<<" "<<Y[i]<<"\n";
	
		}
	fout.close();

}
