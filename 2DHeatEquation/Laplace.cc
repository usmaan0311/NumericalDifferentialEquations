#include"Laplace.h"

Laplace::Laplace()
{}
Laplace::Laplace(int Nx, int Ny, int NST, real xspan, real yspan, real DBC1, real DBC2)
	:Nx(Nx), Ny(Ny), NST(NST), xspan(xspan), yspan(yspan), DBC1(DBC1), DBC2(DBC2)
{

TNX = Nx + 0*NST;//0 Numen BC along x direction
TNY = Ny + 2*NST;//2 Numen BC along y direction

PNX1 = 0;
PNX2 = Nx;

PNY1 = NST;
PNY2 = TNY - NST;

dx = xspan/(Nx - 1);
dy = yspan/(Ny - 1);
std::cout<<"dx"<<dx<<std::endl;
std::cout<<"dy"<<dy<<std::endl;
real IG = 0; 

T = std::vector<std::vector<real>>(TNX, std::vector<real> (TNY, IG ));

x.assign(Nx,0);

for(int i=0; i<Nx;i++)
{
//x.emplace_back(i*dx);
x[i] = i*dx;
}

y.assign(Ny,0);

for(int i=0; i<Ny;i++)
{
//y.emplace_back(i*dy);
y[i] = i*dy;
}


}

real Laplace::l2_norm(std::vector<std::vector<real>> S, std::vector<std::vector<real>> So)
{
real Dsum=0,Usum=0,norm=0;

for(int i=0; i<S.size();i++)
{
for(int j=0; j<S[i].size();j++)
{
	Dsum += pow(So[i][j],2) ;
	Usum += pow( (S[i][j] - So[i][j]),2 );

}

}
norm = sqrt(Usum/Dsum);
//std::cout<<"norm l2\t:"<<norm<<std::endl;
return norm;
}



std::vector<std::vector<real>> Laplace::Jacobi(real tol, int maxiter,std::vector<int>* IterVec, std::vector<real>* ResVec)
{
int iter=0;
real error=1.0;
std::vector<std::vector<real>> To;//(TNX, std::vector<real> (TNY, 0.0 ));
while( (error>tol) && (iter<maxiter) )
{
	To=T;
//	std::cout<<"To size \t:"<<To.size()<<std::endl;
//	std::cout<<"T size \t:"<<T.size()<<std::endl;


	for(int i=PNX1+1;i<PNX2-1;i++)
	{
		for(int j=PNY1;j<PNY2;j++)
		{	
/*	std::cout<<"entering double for 2"<<std::endl;
	std::cout<<"PNX1 :\t"<<PNX1<<std::endl;
	std::cout<<"PNX2 :\t"<<PNX2<<std::endl;
	std::cout<<"PNY1 :\t"<<PNY1<<std::endl;
	std::cout<<"PNY2 :\t"<<PNY2<<std::endl;
	std::cout<<"TNX :\t"<<TNX<<std::endl;
	std::cout<<"TNY :\t"<<TNY<<std::endl;
*/			T[i][j] = ( To[i][j-1]*pow(dx,2) + To[i][j+1]*pow(dx,2) + To[i+1][j]*pow(dy,2) + To[i-1][j]*pow(dy,2)  )/(2*( pow(dx,2) + pow(dy,2) ));
//	std::cout<<"entering after double for 2"<<std::endl;
//			std::cout<<"To["<<i<<"]["<<j<<"] \t :"<<To[i][j]<<std::endl;

		}	

	}	

// Left BC

	for(int j=PNY1;j<PNY2;j++)
	{	
		T[PNX1][j] = 30;
//		std::cout<<" left bc T["<<PNX1<<"]["<<j<<"] \t :"<<T[PNX1][j]<<std::endl;

	}	

// Right BC

	for(int j=PNY1;j<PNY2;j++)
	{
		T[PNX2-1][j] = 60;
//		std::cout<<" right bc T["<<PNX2<<"]["<<j<<"] \t :"<<T[PNX2-1][j]<<std::endl;

	}

// Bottom BC

	for(int i=PNX1;i<PNX2;i++)
	{
		T[i][NST+1] = To[i][NST-1];
//		modified

//		T[i][NST] = ( To[i+1][NST] + To[i-1][NST]  + 2*To[i][NST+1] )/4;
//		std::cout<<" bottom bc T["<<i<<"]["<<NST-1<<"] \t :"<<T[i][NST-1]<<std::endl;

	}

// Top BC

	for(int i=PNX1; i<PNX2;i++)
	{
		T[i][PNY2+NST] = To[i][PNY2-NST] + 2*dy*( To[i][PNY2] - 60 );
//		std::cout<<" top bc T["<<i<<"]["<<PNY2+NST<<"] \t :"<<T[i][PNY2+NST]<<std::endl;

	}

	IterVec->emplace_back(iter);
//	std::cout<<"IterVec["<<0<<"] \t :"<<(*IterVec)[0]<<std::endl;
	ResVec->emplace_back(error);
//	std::cout<<" ResVec["<<0<<"] \t :"<<ResVec->at(0)<<std::endl;
	iter+=1;
	error=l2_norm(T,To);
//	std::cout<<"error l2\t:"<<error<<std::endl;

/*
	for(int i=0; i<T.size();i++)
	{
		for(int j=0; j<T[i].size(); j++)
		{

		std::cout<<"T["<<i<<"]["<<j<<"] \t :"<<T[i][j]<<std::endl;

		}
	}
*/

}


return T;

}

std::vector<std::vector<real>> Laplace::GaussSeidel(real tol, int maxiter,std::vector<int>* IterVec, std::vector<real>* ResVec)
{
int iter=0;
real error=1.0;
std::vector<std::vector<real>> To;//(TNX, std::vector<real> (TNY, 0.0 ));
while( (error>tol) && (iter<maxiter) )
{
	To=T;
//	std::cout<<"To size \t:"<<To.size()<<std::endl;
//	std::cout<<"T size \t:"<<T.size()<<std::endl;


	for(int i=PNX1+1;i<PNX2-1;i++)
	{
		for(int j=PNY1;j<PNY2;j++)
		{	
/*	std::cout<<"entering double for 2"<<std::endl;
	std::cout<<"PNX1 :\t"<<PNX1<<std::endl;
	std::cout<<"PNX2 :\t"<<PNX2<<std::endl;
	std::cout<<"PNY1 :\t"<<PNY1<<std::endl;
	std::cout<<"PNY2 :\t"<<PNY2<<std::endl;
	std::cout<<"TNX :\t"<<TNX<<std::endl;
	std::cout<<"TNY :\t"<<TNY<<std::endl;
*/			T[i][j] = ( T[i][j-1]*pow(dx,2) + T[i][j+1]*pow(dx,2) + T[i+1][j]*pow(dy,2) + T[i-1][j]*pow(dy,2)  )/(2*( pow(dx,2) + pow(dy,2) ));
//	std::cout<<"entering after double for 2"<<std::endl;
//			std::cout<<"To["<<i<<"]["<<j<<"] \t :"<<To[i][j]<<std::endl;

		}	

	}	

// Left BC

	for(int j=PNY1;j<PNY2;j++)
	{	
		T[PNX1][j] = 30;
//		std::cout<<" left bc T["<<PNX1<<"]["<<j<<"] \t :"<<T[PNX1][j]<<std::endl;

	}	

// Right BC

	for(int j=PNY1;j<PNY2;j++)
	{
		T[PNX2-1][j] = 60;
//		std::cout<<" right bc T["<<PNX2<<"]["<<j<<"] \t :"<<T[PNX2-1][j]<<std::endl;

	}

// Bottom BC

	for(int i=PNX1;i<PNX2;i++)
	{
		T[i][NST+1] = T[i][NST-1];
	//	modified
	
	//	T[i][NST] = ( T[i+1][NST] + T[i-1][NST]  + 2*T[i][NST+1] )/4;
//		std::cout<<" bottom bc T["<<i<<"]["<<NST-1<<"] \t :"<<T[i][NST-1]<<std::endl;

	}

// Top BC

	for(int i=PNX1; i<PNX2;i++)
	{
		T[i][PNY2+NST] = T[i][PNY2-NST] + 2*dy*( T[i][PNY2] - 60 );
//		std::cout<<" top bc T["<<i<<"]["<<PNY2+NST<<"] \t :"<<T[i][PNY2+NST]<<std::endl;

	}

	IterVec->emplace_back(iter);
//	std::cout<<"IterVec["<<0<<"] \t :"<<(*IterVec)[0]<<std::endl;
	ResVec->emplace_back(error);
//	std::cout<<" ResVec["<<0<<"] \t :"<<ResVec->at(0)<<std::endl;
	iter+=1;
	error=l2_norm(T,To);
//	std::cout<<"error l2\t:"<<error<<std::endl;

/*
	for(int i=0; i<T.size();i++)
	{
		for(int j=0; j<T[i].size(); j++)
		{

		std::cout<<"T["<<i<<"]["<<j<<"] \t :"<<T[i][j]<<std::endl;

		}
	}
*/

}


return T;

}

std::vector<std::vector<real>> Laplace::SOR(real tol, int maxiter,std::vector<int>* IterVec, std::vector<real>* ResVec, real omega)
{
int iter=0;
real error=1.0;
std::vector<std::vector<real>> To;//(TNX, std::vector<real> (TNY, 0.0 ));
while( (error>tol) && (iter<maxiter) )
{
	To=T;
//	std::cout<<"To size \t:"<<To.size()<<std::endl;
//	std::cout<<"T size \t:"<<T.size()<<std::endl;


	for(int i=PNX1+1;i<PNX2-1;i++)
	{
		for(int j=PNY1;j<PNY2;j++)
		{	
/*	std::cout<<"entering double for 2"<<std::endl;
	std::cout<<"PNX1 :\t"<<PNX1<<std::endl;
	std::cout<<"PNX2 :\t"<<PNX2<<std::endl;
	std::cout<<"PNY1 :\t"<<PNY1<<std::endl;
	std::cout<<"PNY2 :\t"<<PNY2<<std::endl;
	std::cout<<"TNX :\t"<<TNX<<std::endl;
	std::cout<<"TNY :\t"<<TNY<<std::endl;
*/			T[i][j] = (1.0 - omega)*T[i][j] + omega*( T[i][j-1]*pow(dx,2) + T[i][j+1]*pow(dx,2) + T[i+1][j]*pow(dy,2) + T[i-1][j]*pow(dy,2)  )/(2*( pow(dx,2) + pow(dy,2) ));
//	std::cout<<"entering after double for 2"<<std::endl;
//			std::cout<<"To["<<i<<"]["<<j<<"] \t :"<<To[i][j]<<std::endl;

		}	

	}	

// Left BC

	for(int j=PNY1;j<PNY2;j++)
	{	
		T[PNX1][j] = 30;
//		std::cout<<" left bc T["<<PNX1<<"]["<<j<<"] \t :"<<T[PNX1][j]<<std::endl;

	}	

// Right BC

	for(int j=PNY1;j<PNY2;j++)
	{
		T[PNX2-1][j] = 60;
//		std::cout<<" right bc T["<<PNX2<<"]["<<j<<"] \t :"<<T[PNX2-1][j]<<std::endl;

	}

// Bottom BC

	for(int i=PNX1;i<PNX2;i++)
	{
		T[i][NST+1] = T[i][NST-1];
	//	modified
//		T[i][NST] = ( T[i+1][NST] + T[i-1][NST]  + 2*T[i][NST+1] )/4;
//		std::cout<<" bottom bc T["<<i<<"]["<<NST-1<<"] \t :"<<T[i][NST-1]<<std::endl;

	}

// Top BC

	for(int i=PNX1; i<PNX2;i++)
	{
		T[i][PNY2+NST] = T[i][PNY2-NST] + 2*dy*( T[i][PNY2] - 60 );
//		std::cout<<" top bc T["<<i<<"]["<<PNY2+NST<<"] \t :"<<T[i][PNY2+NST]<<std::endl;

	}

	IterVec->emplace_back(iter);
//	std::cout<<"IterVec["<<0<<"] \t :"<<(*IterVec)[0]<<std::endl;
	ResVec->emplace_back(error);
//	std::cout<<" ResVec["<<0<<"] \t :"<<ResVec->at(0)<<std::endl;
	iter+=1;
	error=l2_norm(T,To);
//	std::cout<<"error l2\t:"<<error<<std::endl;

/*
	for(int i=0; i<T.size();i++)
	{
		for(int j=0; j<T[i].size(); j++)
		{

		std::cout<<"T["<<i<<"]["<<j<<"] \t :"<<T[i][j]<<std::endl;

		}
	}
*/

}


return T;

}

void Laplace::PrintSol(std::vector<std::vector<real>> Sol, std::vector<real> Res, std::vector<int> Itr, std::string strg)
{
	if(strg=="Solution")
	{
		for(int i=0; i<Sol.size();i++)
		{
			for(int j=0; j<Sol[i].size(); j++)
			{

				std::cout<<"T["<<i<<"]["<<j<<"] \t :"<<Sol[i][j]<<std::endl;

			}
		}	


	}
	else if(strg=="Iteration")
	{
		for(int i=0; i<Itr.size();i++)
		{
			std::cout<<"Iter["<<i<<"] \t :"<<Itr[i]<<std::endl;

		}
	}
	else if(strg=="Residual")
	{std::cout<<"Inside Residual"<<std::endl;
		for(int i=0; i<Res.size();i++)
		{	
			std::cout<<"Res["<<i<<"] \t :"<<Res[i]<<std::endl;

		}
	}

}

void Laplace::WriteSol(std::vector<std::vector<real>> Sol, std::string fname)
{
	std::ofstream WriteS(fname);
	WriteS.precision(10);
	WriteS.setf(std::ios::scientific);
	assert(WriteS.is_open());
	for(int i=PNX1;i<PNX2;i++)
	{
		for(int j=PNY1; j<PNY2; j++)
		{
		
			WriteS<<Sol[i][j]<<" ";
		}
		WriteS<<"\n";

	}
	WriteS.close();
}
void Laplace::WriteXY(std::string fnamex, std::string fnamey)
{
	std::ofstream WriteX(fnamex);
	WriteX.precision(10);
	WriteX.setf(std::ios::scientific);
	assert(WriteX.is_open());
	for(int i=PNX1;i<PNX2;i++)
	{
		WriteX<<x[i]<<"\n";

	}
	WriteX.close();

	std::ofstream WriteY(fnamey);
	WriteY.precision(10);
	WriteY.setf(std::ios::scientific);
	assert(WriteY.is_open());
	for(int i=PNY1;i<PNY2;i++)
	{
		WriteY<<y[i]<<"\n";

	}
	WriteY.close();
}
void Laplace::WriteResidual(std::vector<real> Res, std::string fname)
{
	std::ofstream WriteRes(fname);
	WriteRes.precision(10);
	WriteRes.setf(std::ios::scientific);
	assert(WriteRes.is_open());
	for(int i=PNX1;i<PNX2;i++)
	{
		WriteRes<<Res[i]<<"\n";

	}
	WriteRes.close();
}
void Laplace::WriteIter(std::vector<int> Itr, std::string fname)
{
	std::ofstream WriteItr(fname);
	WriteItr.precision(10);
	WriteItr.setf(std::ios::scientific);
	assert(WriteItr.is_open());
	for(int i=PNX1;i<PNX2;i++)
	{
		WriteItr<<Itr[i]<<"\n";

	}
	WriteItr.close();
}
