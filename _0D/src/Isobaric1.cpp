#include<iostream>
#include<vector>
#include<fstream>
#include<Eigen/Dense>
#include<cmath>

#include"Jvala_0D.HPP"


double Isobaric::Calc_V(Eigen::VectorXd XFinal) //ch2
{
	double sum=0;	
	for(unsigned i=0;i<N;i++)
		sum+=(XFinal(i)*MW(i));
	double V=m/sum;
	return V;
}

double Isobaric::TempODE(double T,Eigen::VectorXd XFinal,double t)//do sth about t? //ch2
{
	double sum1=0;
	double sum2=0;	
	Eigen::VectorXd Enth,Omega;
	Enth=Eigen::VectorXd::Zero(N);
	Omega=Eigen::VectorXd::Zero(N);
	Enth=Enthalpy(T,t);
	Omega=RateofProduction(T,XFinal,t);
	for(unsigned i=0;i<N;i++)
	{
		sum1+=(Enth(i)*Omega(i));
		sum2+=(XFinal(i)*Calc_Cp(T,i));
	}
	double V=Calc_V(XFinal);
	double value=((Qdot/V)-sum1)/sum2;
	
	return value;
}

double Isobaric::XiODE( double T, Eigen::VectorXd XFinal,double t,int i) //ch1
{
	double sumOmega=0;
	double sumXj=0;
	Eigen::VectorXd Omega;
	Omega=Eigen::VectorXd::Zero(N);
	Omega=RateofProduction(T,XFinal,t);    //chould i store omega instead?
	for(unsigned j=0;j<N;j++)
	{
		sumOmega+=Omega(j);
		sumXj+=XFinal(j);
	}
	double valuei= Omega(i)-XFinal(i)*((sumOmega/sumXj)+(TempODE(T,XFinal,t)/T));
	return valuei;
}
