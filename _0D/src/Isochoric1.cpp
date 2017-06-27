#include<iostream>
#include<vector>
#include<fstream>
#include<Eigen/Dense>
#include<cmath>

#include"Jvala_0D.HPP"




double Isochoric::TempODE( double T, Eigen::VectorXd XFinal, double t)//ch2
{
	double sumOmega=0;
	double sum1=0;
	double sum2=0;	
	Eigen::VectorXd Omega;
	Omega=Eigen::VectorXd::Zero(N);	
	Omega=RateofProduction(T,XFinal,t);		//better to save enth and omega?
	Eigen::VectorXd Enth;
	Enth=Eigen::VectorXd::Zero(N);		
	Enth=Enthalpy(T,t);
	for(unsigned i=0;i<N;i++)
	{
		sumOmega+=Omega(i);
		sum1+=(Enth(i)*Omega(i));
		sum2+=(XFinal(i)*(Calc_Cp(T,i)-R));
	}
	double value=((Qdot/V)+(R*T*sumOmega-sum1))/sum2;
	return value;
}

//double Isochoric::XiODE(Eigen::VectorXd Omega, Eigen::VectorXd XFinal, double T)
//{
//	return Omega(i);
//}											thought it would be better to use Rate of production func or should i save omega(if possible)?

double Isochoric::Calc_P(double T,Eigen::VectorXd XFinal)//ch2
{
	double P=0;	
	for(unsigned i=0;i<N;++i)
		P+=(XFinal(i)*R*T);
	return P;
}

double Isochoric::RatePressure(double T,Eigen::VectorXd XFinal, double t)//ch2
{
	double sumOmega=0;
	double sum2=0;
	Eigen::VectorXd Omega;
	Omega=Eigen::VectorXd::Zero(N);	
	Omega=RateofProduction(T,XFinal,t);	
	Eigen::VectorXd Enth;
	Enth=Eigen::VectorXd::Zero(N);		
	Enth=Enthalpy(T,t);
	for(unsigned i=0;i<N;i++)
	{
		sumOmega+=Omega(i);
		sum2+=(XFinal(i)*TempODE(T,XFinal,t));
	}
	double value=(R*T*sumOmega)+(R*sum2);
	return value;
}
