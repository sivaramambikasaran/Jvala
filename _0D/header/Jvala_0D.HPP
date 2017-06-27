#ifndef __Jvala_0D_HPP__
#define __Jvala_0D_HPP__


#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<Eigen/Dense>
#include<cmath>


//Properties class=========================================================================================================================================================================
class Properties
{
	public:
	//data	
	static constexpr const double R=8.314;
	static constexpr const double Tref=298.15;

	//user input
	double m,Qdot,To; //m-mass of reaction mixture(constant)	Qdot-Rate of heat loss from reactor	To-Initial temperature
	Eigen::VectorXd XInitial;// 	Xinitial-initial conc of species					(should i calc Xinitial from phi?)

	//from files
	unsigned L,N,A; //L-total no. of rxns		N-total no. of species(including other)		A-no. of atoms involved 
	double Tmid,Tupper,Tlower;	//Define temperature ranges for Cp_coefs
	Eigen::MatrixXd reactant_coefs,product_coefs,Cp_coefs; //reactant_coefs-define coefficient matrix for reactants 	product_coefs- for products	Cp_coefs-coefficients to calculate the Cp/H/S of a species at a particular temperature
	Eigen::VectorXd MW,Af,Eaf,nf,Ar,Ear,nr; //MW-Molecular weights	 Af-preexponential factor for forward rxns	Eaf-activation energ for forward rxns	nr-Exponent on temperature in the expression to calcute rate constant		Ar,Ear,nr-for reverse rxns
	std::vector<std::string> Species,Atoms; // species and atoms
	std::vector<std::vector<double> > Press_coefs; //Coefficients to dtermine rate constant for pressure dependent rxn rates.
	Eigen::MatrixXd Mxx; //Matrix to store the composition of third bodies Ms
	std::vector<std::vector<unsigned> > M_index; //stores the indices of rxns with third bodies

/*	to calc
	std::vector<Eigen::VectorXd> XFinalAll;
	Eigen::VectorXd TAll;
	std::vector<Eigen::VectorXd> KfAll,KrAll;
	std::vector<Eigen::VectorXd> EFinalAll;//all enthalpies
	std::vector<Eigen::VectorXd> OmegaAll;
*/	
	
	//CONSTRUCTOR
	Properties(double To, Eigen::VectorXd Xi,double Q, double m);	//works

	//functions 
	void resize_vectors();						//ok
	double Calc_K(double A, double n, double E,double T,Eigen::VectorXd Xi, unsigned i);														//works
	double Calc_H(double T, unsigned i);				//works
	double Calc_S(double T, unsigned i);				//works
	double Calc_Cp(double T,unsigned i);				//works
	double Calc_Kb(double T,Eigen::VectorXd Xi,unsigned k);		//works
	double Calc_F(double P_r,double T,unsigned i);			//works
	double Calc_M(Eigen::VectorXd Xi,unsigned i);			//works
	void Calc_MW();	//to get MW					//ok
	int get_Index(std::string s,std::vector<std::string> Sp);	//ok
	void get_Cpcoefs(); //to get Cp_coefs				//ok
	void get_Species();						//ok
	void get_rpcoefs();//one to get a,n,E,coef matrices vectors from .mech file for use in Calc_K //ok
	void thirdbodies();						
	Eigen::VectorXd Calc_RateProgress(double T,Eigen::VectorXd XFinal, double t);										
	Eigen::VectorXd RateofProduction(double T,Eigen::VectorXd XFinal,double t);//to calc rate at a particular t							
	Eigen::VectorXd Enthalpy(double T, double t);//at a particular T												
};

//Isobaric class=========================================================================================================================================================================
class Isobaric : public Properties
{
	// do i keep an array for all calculated V?	
	public:
	Isobaric(double To, Eigen::VectorXd Xi,double Q, double m) :Properties(To,Xi,Q,m)										
	{
		
	};
	double Calc_V(Eigen::VectorXd XFinal);		//returns V/mol fuel												
	double TempODE(double T, Eigen::VectorXd XFinal,double t);													
	double XiODE(double T, Eigen::VectorXd XFinal,double t,int i);						
};
//Isochoric class=========================================================================================================================================================================
class Isochoric : public Properties
{
	public:
	double V;
	Isochoric(double To, Eigen::VectorXd Xi,double Q, double m, double V):Properties(To,Xi,Q,m)		
	{
		this->V=V;
	};
	double Calc_P(double T,Eigen::VectorXd XFinal);							
	double RatePressure(double T,Eigen::VectorXd XFinal, double t);						
	double TempODE( double T, Eigen::VectorXd XFinal, double t);						
//	double XiODE(Eigen::VectorXd Omega, Eigen::VectorXd XFinal, double T);
};
#endif //(__Jvala_0D_HPP__)
