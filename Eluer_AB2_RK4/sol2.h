#ifndef sol2
#define sol2
#include<vector>
#include<string>

typedef double real;

class Sol2{
real ic;
real h;
real x0;
real xf;
int n;
std::vector<real> T,Y;


	public:
		Sol2();
		Sol2(real , real, real, real);
		std::vector<std::string>	ReadParams(std::string);
		std::vector<real>	Euler(real (*)(real,real));
		std::vector<real>	AB2(real (*)(real,real));
		std::vector<real>	RK4(real (*)(real,real));
		std::vector<std::vector<real>>	AvgErr(real (*)(real), real (*)(real, real), int);
		void	WriteSol(std::vector<real>&, std::string, int);
		void	WriteErr(std::vector<real>&, std::string, int);

};

#endif
