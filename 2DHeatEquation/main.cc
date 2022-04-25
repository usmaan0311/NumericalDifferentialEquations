#include"Laplace.cc"

int main(int argc, char** argv)
{
int Nx=128, Ny=256, NST=1, maxiter=20000000;
real xspan=1, yspan=2, DBC1=30, DBC2=60, error=1e-6, omega= ( 2.0/(1.0 + ( M_PI/( std::max(Nx,Ny) ) ) ) );

Laplace L;
Laplace L1(Nx,Ny,NST,xspan,yspan,DBC1,DBC2);

std::vector<std::vector<real>> SV;
std::vector<real> Res;
std::vector<int> Itr;

auto ti = std::chrono::high_resolution_clock::now();

//SV=L1.Jacobi(error, maxiter, &Itr, &Res);
SV=L1.GaussSeidel(error, maxiter, &Itr, &Res);
//SV=L1.SOR(error, maxiter, &Itr, &Res, omega);

auto tf = std::chrono::high_resolution_clock::now();

auto dt = std::chrono::duration_cast<std::chrono::microseconds>(tf - ti);

std::thread t0(&Laplace::PrintSol,&L1,SV,Res,Itr,"Solution");
std::thread t01(&Laplace::PrintSol,&L1,SV,Res,Itr,"Residual");



std::thread t1(&Laplace::WriteSol,&L1,SV,"Solution.txt");
std::thread t2(&Laplace::WriteXY,&L1,"X.txt","Y.txt");
std::thread t3(&Laplace::WriteResidual,&L1,Res,"Residual.txt");
std::thread t4(&Laplace::WriteIter,&L1,Itr,"Iteration.txt");
t0.join();
t01.join();
t1.join();
t2.join();
t3.join();
t4.join();

std::cout<<"Solution done in \t: "<<dt.count()/1000<<" ms"<<std::endl;
std::cout<<"Solution converged in \t: "<<Itr.back()<<" iterations"<<std::endl;
std::cout<<"Final residual \t: "<<Res.back()<<std::endl;

//L1.WriteSol(SV,"Solution.txt");
return EXIT_SUCCESS;
}
