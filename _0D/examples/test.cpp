#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include<Eigen/Dense>

#include"anukalana.hpp"
#include"Jvala_0D.HPP"

double tInit				=	0.0;
double tFinal				=	1000;
long unsigned nTimeSteps;
double deltat;
int No	=	35;



class myIntegrator: public Integrator
{
	protected:
	Isobaric Problem;

	Eigen::VectorXd function(double t, Eigen::VectorXd y)
	{
		Eigen::VectorXd temp(y.rows());
		for(unsigned i=0;i<N;i++)
		{
			temp(i)=Problem.XiODE(y(N),y,t,i);
		}
		temp(N)	=Problem.TempODE(y(N),y,t); 
		
		return temp;
	};


	public:
	myIntegrator(Eigen::VectorXd yInitial, double tInit, int nTimeSteps, double deltat,Isobaric Prob) : Integrator(yInitial, tInit, nTimeSteps, deltat)
	{
		this->Problem=Prob;
	};
	~myIntegrator(){};
};

void writeToFile(std::string fileName, std::vector<Eigen::VectorXd> yFinalComputedAll) 
{
	std::ofstream myfile;
	myfile.open(fileName.c_str(),std::ios::out);
	myfile << "\\addplot[mark = o,red] coordinates {";
	for (unsigned long int j=0; j<nTimeSteps; j+=nTimeSteps/1000) 
	{
		myfile << "(" << tInit+j*deltat ;
		for(unsigned k=0;k<No;++k)
			myfile  << ","<< yFinalComputedAll[j](k);
		myfile << ")";
	}
	myfile << "};\n";
	// myfile << "\\addplot[mark = o,blue] coordinates {";
	// for (unsigned long int j=0; j<nTimeSteps; j+=nTimeSteps/1000) {
	// 	myfile << "(" << tInit+j*deltat << "," << yFinalComputedAll[j](1) << ")";
	// }
	// myfile << "};\n";
	myfile.close();
}

void writeToFile(std::string fileName, std::vector<double>& timesteps, std::vector<Eigen::VectorXd> yFinalComputedAll) 
{
	std::ofstream myfile;
	myfile.open(fileName.c_str(),std::ios::out);
	myfile << "\\addplot[mark = o,red] coordinates {";
	long unsigned int nTimeSteps	=	timesteps.size();
	for (unsigned long int j=0; j<nTimeSteps; j+=nTimeSteps/1000) 
	{
		myfile << "(" << timesteps[j];
		for(unsigned k=0;k<No;++k)
			myfile  << ","<< yFinalComputedAll[j](k);
		myfile << ")";
	}
	myfile << "};\n";
	// myfile << "\\addplot[mark = o,blue] coordinates {";
	// for (unsigned long int j=0; j<nTimeSteps; j+=nTimeSteps/1000) {
	// 	myfile << "(" << timesteps[j] << "," << yFinalComputedAll[j](1) << ")";
	// }
	// myfile << "};\n";
	myfile.close();
}

int main(int argc, char* argv[]) 
{
	nTimeSteps	=	double(atoi(argv[1]));
	deltat		=	double(tFinal-tInit)/double(nTimeSteps);
	srand(time(NULL));

	double To=500.0; //random
	Eigen::VectorXd Xi(34);
	for(unsigned i=0;i<34;++i) Xi(i)=1.0; //random
	double Q=0.0;//random
	double m=0.0;	//random

	Isobaric P;
	std::string temp;
	Eigen::VectorXd yInitial(P.N);
	yInitial=Eigen::VectorXd::Zero(P.N);
	std::ifstream fin;
	fin.open("./inputFiles/allInput.txt",std::ios::in);
	for(unsigned i=0;i<10;++i)
		std::getline(fin,temp);	

	while(fin)
	{
		fin>>temp;
		int i=P.get_Index(temp,P.Species);
		fin>>temp;
		if(i==-1) continue;
		yInitial(i)=std::atof(temp.c_str());
	}
	myIntegrator(yInitial,tInit,nTimeSteps,deltat,P);
	

/*
	Eigen::VectorXd Enth,Omega;
	Enth=Eigen::VectorXd::Zero(A.N);
	Omega=Eigen::VectorXd::Zero(A.N);
	Enth=A.Enthalpy(To,3);
	Omega=A.RateofProduction(To,Xi,3);
	for(unsigned j=0;j<A.N;++j)
	{
		std::cout<<A.Species[j]<<"\nEn= "<<Enth(j)<<"\tOm= "<<Omega(j)<<"\n";
	}
*/
	//air=80N2 and 20O2. Fuel=
/*
	unsigned i=6;
	for(unsigned j=0;j<A.N;++j)
	{
		if(A.reactant_coefs(i,j)!=0)
		std::cout<<A.Species[j]<<" + ";
	}
	std::cout<<"-->";
	for(unsigned j=0;j<A.N;++j)
	{
		if(A.product_coefs(i,j)!=0)
		std::cout<<A.Species[j]<<" + ";
	}
	std::cout<<"rate coesf: "<<A.Af(i)<<"\t"<<A.nf(i)<<"\t"<<A.Eaf(i)<<"\n"<<A.Ar(i)<<"\t"<<A.nr(i)<<"\t"<<A.Ear(i)<<"\n";
	std::cout<<"Pressure coefs: ";
	unsigned t;
	for(t=0;t<A.Press_coefs.size();++t)
	{
		if(i==A.Press_coefs[t][0]) break;
	}
	for(unsigned j=0;j<9;++j)
		std::cout<<A.Press_coefs[t][j]<<"\t";
	std::cout<<"\nCp coefs & MW: \n";
	
	for(unsigned j=0;j<A.N;++j)
	{
		if(A.reactant_coefs(i,j)!=0)
		{
		std::cout<<A.Species[j]<<"\t"<<A.MW(j)<<"\n";
		for(unsigned k=0;k<14;++k)
		std::cout<<k<<":"<<A.Cp_coefs(j,k)<<"\t";
		}
	}
	for(unsigned j=0;j<A.N;++j)
	{
		if(A.product_coefs(i,j)!=0)
		{
		std::cout<<A.Species[j]<<"\t"<<A.MW(j)<<"\n";
		for(unsigned k=0;k<14;++k)
		std::cout<<k<<":"<<A.Cp_coefs(j,k)<<"\t";
		}
	}
*/

	std::vector<Eigen::VectorXd> yFinalComputedAll;

	std::ofstream myfile;
	std::string fileName;

	//	EulerExplicit or RK1
	std::vector<double> parameters;
	yFinalComputedAll	=	A.RK_NonAdaptive_All(1,1,parameters);
	fileName			=	"./Output/Isobaric_EulerExplicit_NonAdaptive_All.txt";
	writeToFile(fileName, yFinalComputedAll);
/*
	//	RK2
	parameters.clear();
	parameters.push_back(0.5);
	yFinalComputedAll	=	A.RK_NonAdaptive_All(2,1,parameters);
	fileName			=	"./Output/Isobaric_RK2stage_NonAdaptive_All.txt";
	writeToFile(fileName, yFinalComputedAll);

*/
	
	return 0;
	
}
