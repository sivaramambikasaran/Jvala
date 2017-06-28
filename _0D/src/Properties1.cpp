#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<Eigen/Dense>
#include<cmath>
#include<cctype>
#include<regex>

#include"Jvala_0D.HPP"

std::string mechanismfile,thermofile,transportfile;


void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

double Properties::Calc_K(double A, double n, double E,double T,Eigen::VectorXd Xi, unsigned i)
{
	//units?!
std::cout<<T<<"\n";
	if(A==0&&n==0&&E==0)
		return Calc_Kb(T,Xi,i);
	double K_inf= A*(std::pow(T,n))*(std::exp(-(E*1000.0)/(R*T)));
	unsigned j;
	for(j=0;j<Press_coefs.size();++j)
	{
		if(Press_coefs[j][0]==i) break;
	}
	if(j==Press_coefs.size())
		return K_inf;

	double K_o= Press_coefs[j][1]*(std::pow(T,Press_coefs[j][2]))*(std::exp(-(Press_coefs[j][3]*1000.0)/(R*T)));
	double P_r=(K_o/K_inf)*Calc_M(Xi,i);
	double K=Calc_F(P_r,T,j)*(K_inf/(1+((K_inf/K_o)*Calc_M(Xi,i))));
	return K;
}

double Properties::Calc_M(Eigen::VectorXd Xi,unsigned i)
{
	//give a marker to rxns with third bodies.
	//if thrid bodies there, compute [M] else sum [Xi] of the reaction
	unsigned k;
	for(k=0;k<M_index.size();++k)
	{
		for(unsigned j=0;j<M_index[k].size();++j)
		{
			if(M_index[k][j]==i)
				break;
		}
	}
	if(k==M_index.size())
	{
		double sum=0;
		for(unsigned l=0;l<N;++l)
			sum+=Xi(l);
		return sum;
	}
	else
	{
		double sum=0;
		for(unsigned l=0;l<N;++l)
			sum+=Mxx(i,l)*Xi(i);
		return sum;
	}

	return 0;
}

double Properties::Calc_H(double T,unsigned i)
{
	double H=9;
	if((T>Tmid)&&(T<=Tupper))
	H=Cp_coefs(i,0) + Cp_coefs(i,1)*(T/2) + Cp_coefs(i,2)*std::pow(T,(2/3)) + Cp_coefs(i,3)*std::pow(T,(3/4)) + Cp_coefs(i,4)*std::pow(T,(4/5)) + Cp_coefs(i,5)/T;
	else if((T>=Tlower)&&(T<=Tmid))
	H=Cp_coefs(i,7) + Cp_coefs(i,8)*(T/2) + Cp_coefs(i,9)*std::pow(T,(2/3)) + Cp_coefs(i,10)*std::pow(T,(3/4)) + Cp_coefs(i,11)*std::pow(T,(4/5)) + Cp_coefs(i,12)/T;
	return T*R*H;
}

double Properties::Calc_S(double T,unsigned i)
{
	double S=6;
	if((T>Tmid)&&(T<=Tupper))
	S=Cp_coefs(i,0)*std::log(T) + Cp_coefs(i,1)*(T) + Cp_coefs(i,2)*std::pow(T,(2))*(1/2) + Cp_coefs(i,3)*std::pow(T,(3))*(1/3) + Cp_coefs(i,4)*std::pow(T,(4))*(1/4) + Cp_coefs(i,6);
	else if((T>=Tlower)&&(T<=Tmid))
	S=Cp_coefs(i,7)*std::log(T) + Cp_coefs(i,8)*(T) + Cp_coefs(i,9)*std::pow(T,(2))*(1/2) + Cp_coefs(i,10)*std::pow(T,(3))*(1/3) + Cp_coefs(i,11)*std::pow(T,(4))*(1/4) + Cp_coefs(i,13);
	return R*S;
}

double Properties::Calc_Cp(double T,unsigned i)
{
	double Cp=8;
	if((T>Tmid)&&(T<=Tupper))
	Cp=Cp_coefs(i,0) + Cp_coefs(i,1)*(T) + Cp_coefs(i,2)*std::pow(T,2) + Cp_coefs(i,3)*std::pow(T,(3)) + Cp_coefs(i,4)*std::pow(T,(4));
	else if((T>=Tlower)&&(T<=Tmid))
	Cp=Cp_coefs(i,7) + Cp_coefs(i,8)*(T) + Cp_coefs(i,9)*std::pow(T,(2)) + Cp_coefs(i,10)*std::pow(T,(3)) + Cp_coefs(i,11)*std::pow(T,(4));
	return R*Cp;
}

double Properties::Calc_Kb(double T,Eigen::VectorXd Xi,unsigned k)
{
	double deltaHr=0;
	double deltaSr=0;
	for(unsigned i=0;i<N;++i)
	{
		if(reactant_coefs(k,i)!=0)
		{
			deltaHr-=(reactant_coefs(k,i)*Calc_H(T,i)); //Tref or T?
			deltaSr-=(reactant_coefs(k,i)*Calc_S(T,i));
		}
		if(product_coefs(k,i)!=0)
		{
			deltaHr+=(product_coefs(k,i)*Calc_H(T,i));
			deltaSr+=(product_coefs(k,i)*Calc_S(T,i));
		}
	}
	double deltaG=deltaHr-deltaSr*T;
	double Keqm=std::exp((-deltaG/(R*T)));
	double Kb=Calc_K(Af(k),nf(k),Eaf(k),T,Xi,k)/Keqm;
	return Kb;
}

double Properties::Calc_F(double P_r,double T,unsigned i)
{
	double F_cent=( Press_coefs[i][4])*std::exp(-T/ Press_coefs[i][5]) +Press_coefs[i][6]*std::exp(-T/Press_coefs[i][7]) + std::exp(-Press_coefs[i][9]/T);
	double C=-0.4-0.67*std::log(F_cent);
	double N=0.75-1.27*std::log(F_cent);
	
	double denom=1+std::pow(((std::log(P_r)+C)/(N-0.14*(std::log(P_r)+C))),2);
	double F= exp(std::log(F_cent)/denom);
	return F;
}

void Properties::Calc_MW()
{
std::cout<<"in MW\n";
	//to find each atoms MW
	std::ifstream fin;
	fin.open("./inputFiles/MolarMasses.txt", std::ios::in);
	Eigen::VectorXd AMW(A);
	AMW=Eigen::VectorXd::Zero(A);
	for(unsigned i=0;i<A;++i)
	{
		std::string token;
		std::string s=Atoms[i];
		fin.seekg(0,std::ios::beg);
		while(fin)
		{
			fin>>token;
			fin>>token;
			for (auto & c: token) c = toupper(c);
			if(token==s)
			{
				fin>>token;
				fin>>token;
				AMW(i)=std::atof(token.c_str());
				break;
			}
		}
std::cout<<AMW(i)<<"\t";
	}
std::cout<<"^atom's MM\n";
	fin.close();
/*Te different molecules
NC7KET35-C7H14O3 - first part not counted
C7H14OOH3-5O2-C7H15O4 ?
C7H14OOH3-5-C7H15O2
N-C7H16 N not counted
P-C4H9

*/
	//to assign MW to each Spc
	std::string ch1="";
	std::string ch2="";
	std::string num="";
	int n=0;
	unsigned j=0;
	MW.resize(N);
	MW=Eigen::VectorXd::Zero(N);
	for(unsigned i=0;i<Species.size();++i)
	{
std::cout<<"k1\t";
		ch1="";
		ch2="";
		num="";
		n=0;
		j=0;

std::cout<<i<<"\t"<<Species[i]<<"\t";
		if(Species[i]=="OTHER")
			continue; //get definition of other
		if(Species[i].find_last_of("-")!=std::string::npos)
			j=Species[i].find_last_of("-")+1;
		while(j<Species[i].length()) 
		{
			if(isalpha(Species[i][j]))
			{
				num="";
				ch1=ch2=Species[i][j];
				if(j<Species[i].length()-1)
					ch2+=Species[i][++j];
				int c1=get_Index(ch1,Atoms);
				int c2=get_Index(ch2,Atoms);
				if(c1!=-1||c2!=-1)
				{
					if(c1!=-1&&c2==-1)
						c2=c1;
					else
						j++;
					
					while(j<Species[i].length()&&isdigit(Species[i][j]))
						{num+=Species[i][j];j++;}
					if(num=="")
						n=1;
					else
						n=std::atof(num.c_str());
					MW(i)+=(n*AMW(c2));
				}
			}
			else
				j++;
		}
std::cout<<MW(i)<<"\t";
std::cout<<"k2\n";
	}
std::cout<<"\nMWs\n";

}

void Properties::get_Cpcoefs()
{
	std::ifstream fin;
	fin.open(thermofile.c_str(),std::ios::in);

	std::string temp;
	std::getline(fin,temp);
	fin>>Tlower>>Tmid>>Tupper;
	std::getline(fin,temp);
	std::size_t start=fin.tellg();

std::cout<<"in get Cp coefs\n";
	int dash=0;
	std::regex Spc_ignore("[A-Z]\\-[A-Z][A-Za-z0-9\\-]*");
	for(unsigned i=0;i<N;++i)
	{
		fin.seekg(start,std::ios::beg);
		std::string s=Species[i];
		dash=0;
		for(unsigned j=0;j<s.length();++j)
		{
			if(s[j]=='-')
				dash++;
		}
		if(dash==2)
			s=s.std::string::substr(0,s.std::string::find_last_of("-"));
		else if(dash==1)
		{
			if(!std::regex_match(s,Spc_ignore))
				s=s.std::string::substr(0,s.std::string::find_last_of("-"));	
		}
std::cout<<s<<"\n";

		while(fin)
		{
			std::getline(fin,temp,' ');
			if(temp!=s)
			{
				std::getline(fin,temp);	
				std::getline(fin,temp);	
				std::getline(fin,temp);	
				std::getline(fin,temp);	
				continue;
			}
			else if(temp==s)
			{
				std::getline(fin,temp);
				std::getline(fin,temp);
				Cp_coefs(i,0)=std::atof(temp.std::string::substr(0,15).c_str());
				Cp_coefs(i,1)=std::atof(temp.std::string::substr(15,15).c_str());
				Cp_coefs(i,2)=std::atof(temp.std::string::substr(30,15).c_str());
				Cp_coefs(i,3)=std::atof(temp.std::string::substr(45,15).c_str());
				Cp_coefs(i,4)=std::atof(temp.std::string::substr(60,15).c_str());
	
				std::getline(fin,temp);
				Cp_coefs(i,5)=std::atof(temp.std::string::substr(0,15).c_str());
				Cp_coefs(i,6)=std::atof(temp.std::string::substr(15,15).c_str());
				Cp_coefs(i,7)=std::atof(temp.std::string::substr(30,15).c_str());
				Cp_coefs(i,8)=std::atof(temp.std::string::substr(45,15).c_str());
				Cp_coefs(i,9)=std::atof(temp.std::string::substr(60,15).c_str());
	
				std::getline(fin,temp);
				Cp_coefs(i,10)=std::atof(temp.std::string::substr(0,15).c_str());
				Cp_coefs(i,11)=std::atof(temp.std::string::substr(15,15).c_str());
				Cp_coefs(i,12)=std::atof(temp.std::string::substr(30,15).c_str());
				Cp_coefs(i,13)=std::atof(temp.std::string::substr(45,15).c_str());
				break;
			}
		}
for(unsigned j=0;j<14;++j)
std::cout<<Cp_coefs(i,j)<<"\t";
std::cout<<"\n";
	}
	fin.close();
}

int Properties::get_Index(std::string s,std::vector<std::string> Sp)
{
	for(unsigned i=0;i<Sp.size();++i)
	{
		if(s.compare(Sp[i])==0)
			return i;
	}
	return -1;
}

Properties::Properties() //somewhat a constructor :/
{
	std::ifstream fin;
	fin.open("./inputFiles/allInput.txt",std::ios::in);
	std::string temp;
	std::getline(fin,temp);
	std::getline(fin,mechanismfile);
	std::getline(fin,temp);
	std::getline(fin,temp);
	std::getline(fin,thermofile);
	std::getline(fin,temp);
	std::getline(fin,temp);
	std::getline(fin,transportfile);f
/*
	this->To=To;
	this->XInitial=Xi;
	this->Qdot=Q;
	this->m=m;	
	this->L=1;
	this->N=2;
	this->A=0;
	this->reactant_coefs=Eigen::MatrixXd::Zero(L,N);
	this->product_coefs=Eigen::MatrixXd::Zero(L,N);
	this->Cp_coefs=Eigen::MatrixXd::Zero(N,14);
	this->MW=Eigen::VectorXd::Zero(N);
	this->Af=Eigen::VectorXd::Zero(L);
	this->Eaf=Eigen::VectorXd::Zero(L);
	this->nf=Eigen::VectorXd::Zero(L);
	this->Ar=Eigen::VectorXd::Zero(L);
	this->Ear=Eigen::VectorXd::Zero(L);
	this->nr=Eigen::VectorXd::Zero(L);
*/
	get_Species(); //gets L,N,A,Species
	resize_vectors();
	get_rpcoefs(); //reactant_coefs,product_coefs,A,n,Ea,Press_coefs,M_index,Mxx
	get_Cpcoefs();
	Calc_MW();
}

void Properties::resize_vectors() //call in the input from file function
{
std::cout<<"in resize vectors\n";
	this->reactant_coefs.resize(L,N);
	this->product_coefs.resize(L,N);
	this->Cp_coefs.resize(N,14);
	this->MW.resize(N);
	this->Af.resize(L);
	this->Eaf.resize(L);
	this->nf.resize(L);
	this->Ar.resize(L);
	this->Ear.resize(L);
	this->nr.resize(L);
	this->reactant_coefs=Eigen::MatrixXd::Zero(L,N);
	this->product_coefs=Eigen::MatrixXd::Zero(L,N);
	this->Cp_coefs=Eigen::MatrixXd::Zero(N,14);
	this->MW=Eigen::VectorXd::Zero(N);
	this->Af=Eigen::VectorXd::Zero(L);
	this->Eaf=Eigen::VectorXd::Zero(L);
	this->nf=Eigen::VectorXd::Zero(L);
	this->Ar=Eigen::VectorXd::Zero(L);
	this->Ear=Eigen::VectorXd::Zero(L);
	this->nr=Eigen::VectorXd::Zero(L);


}

Eigen::VectorXd Properties::Calc_RateProgress(double T,Eigen::VectorXd XFinal, double t)//ch2
{	
	Eigen::VectorXd q;
	Eigen::VectorXd Kfor,Krev;
	q=Eigen::VectorXd::Zero(L);
	Kfor=Eigen::VectorXd::Zero(L);
	Krev=Eigen::VectorXd::Zero(L);
	for(unsigned i=0;i<L;i++)
	{
		double XiProduct1=1;
		double XiProduct2=1;
		Kfor(i)=Calc_K(Af(i),nf(i),Eaf(i),T,XFinal,i);
		Krev(i)=Calc_K(Ar(i),nr(i),Ear(i),T,XFinal,i);		//not necessary to create vector if i'm not storing Kf,Kr
		q(i)=0;
		for(unsigned j=0;j<N;j++)
		{
			XiProduct1*=std::pow(XFinal(j),reactant_coefs(i,j));
			XiProduct2*=std::pow(XFinal(j),product_coefs(i,j));
		}
		q(i)=(Kfor(i)*XiProduct1)-(Krev(i)*XiProduct2);
	}
//	KfAll.pushback(Kfor);
//	KrAll.pushback(Krev);
	return q;
	//Check if youre to return reference or by value
}

Eigen::VectorXd Properties::RateofProduction(double T,Eigen::VectorXd XFinal,double t) //ch2
{
	Eigen::VectorXd Omega,q;
	Omega=Eigen::VectorXd::Zero(N);
	q=Eigen::VectorXd::Zero(L);
	q=Calc_RateProgress(T,XFinal,t);
	for(unsigned j=0;j<N;j++)
	{
		Omega(j)=0;
		for(unsigned i=0;i<L;i++)
		{
			Omega(j)+=((product_coefs(i,j)-reactant_coefs(i,j))*q(i));
		}
	}
	return Omega;
}

Eigen::VectorXd Properties::Enthalpy(double T, double t)//what about t? //ch2
{
	Eigen::VectorXd E;
	E=Eigen::VectorXd::Zero(N);	
	for(unsigned i=0;i<N;i++)
//		getindex
		E(i)=Calc_H(T,i);
	return E;
}


void Properties::get_Species()
{
std::cout<<"in get_species\n";

	L=0;
	N=0;
	A=0;

	std::ifstream fin;
	fin.open(mechanismfile.c_str(),std::ios::in);
	std::regex Spc("[A-Z][A-Za-z0-9\\-]*");
	std::regex coef("[0-9]\\.*");
	std::regex line("[L|l]et");
	std::regex Mspc("\\[[A-Za-z][A-Za-z0-9\\-]*\\]");
	std::size_t start,end ;
	std::string token;
if(!fin)
std::perror("Erorr");		
	while(fin)
	{

		std::string s;
		std::getline(fin,s);
		if(s.std::string::length()==0) continue;

		start = 0;
    		end = s.std::string::find_first_of(" \t\n",start);
		token = s.substr(start, end - start);
		if(std::regex_match(token,line))
		{
			start = end + 1;
			end = s.std::string::find_first_of(" \t\n",start);
			std::string token = s.substr(start, end - start);
			if(token=="allowed")
			//get allowed atoms 
			{

				start=s.std::string::find("be");
				end = s.std::string::find_first_of(" \t\n",start);
				do
				{
					start=end+1;
					while(s[start]==' ')
					{
							start=start+1;
					}
					end = s.std::string::find_first_of(" \t\n",start);
					if(end==std::string::npos)
						end=s.length();
					std::string token = s.substr(start, end - start-1);
					for (auto & c: token) c = toupper(c);
std::cout<<token<<"\n";
					Atoms.push_back(token);
					A++;

				}while(end<s.length());
			}
			else if(token=="additional")
			{
std::cout<<"in addi\n";
				start=s.std::string::find("be");
				end = s.std::string::find_first_of(" \t\n",start);
				std::string token = s.substr(start, end - start-1);
std::cout<<token<<"\n";
				do
				{
					start=end+1;
					while(s[start]==' ')
					{
							start=start+1;
					}
					end = s.std::string::find_first_of(" \t\n",start);
					if(end==std::string::npos)
						end=s.length();
					std::string token = s.substr(start, end - start-1);
std::cout<<token<<"\n";
					Species.push_back(token);
					N++;

				}while(end<s.length());
			}

			else if(token[0]=='M')
			{
std::cout<<"in M\n";
				int lch=0;
				while (lch!=1)
    				{
					start = end + 1;
					end = s.std::string::find_first_of(" \t\n",start);
					if(end==std::string::npos)
					{
						lch=1;
						end=s.length()-1;
					}
					std::string token = s.substr(start, end - start);
std::cout<<token<<"\t";
					if(std::regex_match(token,Mspc))
					{
						token=token.std::string::substr(1,end-start-2);
std::cout<<token<<"\n";
						if(get_Index(token,Species)==-1||N==0)
						{
							Species.push_back(token);
							N++;
						}
					}
					else if(std::regex_match(token,coef)||token=="="||token=="+") 
						continue;
    				}
				
			}
			continue;
		}

		if(isdigit(token[0])||token[0]=='C')
		{
			unsigned i=0;
			if(!isdigit(token[0])) i++;
			std::string temp=token.substr(i,1);
			i++;
			for(;i<token.length();++i)
			{
				if(isdigit(token[i])) 
					{temp+=token[i];}
			}
			L=std::atof(temp.c_str());
		}
		else 
			continue;

		while (end != std::string::npos)
    		{
			start = end +1;
			end = s.std::string::find_first_of(" \t\n",start);
			std::string token = s.substr(start, end - start);
			if(std::regex_match(token,Spc))
			{ 
				if(get_Index(token,Species)==-1||N==0)
				{
					Species.push_back(token);
					N++;
				}
			}
			else if(std::regex_match(token,coef)||token=="->"||token=="+") 
				{continue;}
			else if(token=="{")
				{break;}
    		}
	}
std::cout<<"Species size"<<Species.size()<<"\n";
	std::cout<<"L="<<L<<"\nN="<<N<<"\n"<<"A="<<A<<"\n";
	for(int l=0;l<N;++l)
	std::cout<<Species[l]<<"\t";
	std::cout<<"\n";
	for(int l=0;l<A;++l)
	std::cout<<Atoms[l]<<"\t";
	std::cout<<"\n";
	fin.close();
}

void Properties::get_rpcoefs()
{
std::cout<<"in get rpcoefs\n";
	std::ifstream fin;
	fin.open(mechanismfile.c_str(),std::ios::in);
	std::regex Spc("[A-Z][A-Za-z0-9\\-]*");
	std::regex coef("[0-9]*");
	std::regex Mspc("\\[[A-Z][A-Za-z0-9\\-]*\\]");
	std::regex Mline("[L|l]et M%");
	std::size_t start,end ;
	std::string token;

//to get rxn num and f or b
	while(fin)
	{
		std::string s;
		std::getline(fin,s);
		if(s.std::string::length()==0) continue;

		start = 0;
		std::size_t old=0;
    		end = s.std::string::find_first_of(" \t\n",start);
		token = s.substr(start, end - start);
		int k=0;
		char check=' ';

		if(isdigit(token[0])||token[0]=='C')
		{
			unsigned i=0;
			if(!isdigit(token[0])) i++;
			std::string temp=token.substr(i,1);
			i++;
			for(;i<token.length();++i)
			{
				if(isdigit(token[i])) 
					{temp+=token[i];}
				else 
					break;
			}
			check=token[i];
			k=std::atof(temp.c_str());
		}
		else
			continue;

//to get reactant coef matrix
//for f rxns only
		if(check!='b')
		{
   			while (end != std::string::npos)
    			{
				old=start;
				start = end + 1;
				end = s.std::string::find_first_of(" \t\n",start);
				token = s.substr(start, end - start);
				if(token=="->") break;
				int Spi;
				Spi=get_Index(token,Species);
				if(std::regex_match(token,Spc))
				{
					std::string prev=s.std::string::substr(old,start-old-1);
					if(std::regex_match(prev,coef))
					{
						reactant_coefs(k-1,Spi)=std::atof(prev.c_str());
					}
					else
						{reactant_coefs(k-1,Spi)=1;}
				}
				else if(std::regex_match(token,coef)||token=="+") 
					continue;
    			}
std::cout<<"for:"<<k<<check<<"\n";
for(int l=0;l<N;l++)
std::cout<<reactant_coefs(k-1,l)<<"\t";
std::cout<<"\n";

//to get product coef matrix
			while (end != std::string::npos)
	    		{
				old=start;
				start = end + 1;
				end = s.std::string::find_first_of(" \t\n",start);
				token = s.substr(start, end - start);
				if(token=="{") break;
				int Spi;
				Spi=get_Index(token,Species);
				if(std::regex_match(token,Spc))
				{
					std::string prev=s.std::string::substr(old,start-old-1);
					if(std::regex_match(prev,coef))
					{
						product_coefs(k-1,Spi)=std::atof(prev.c_str());
					}
					else
						product_coefs(k-1,Spi)=1;
				}
				else if(std::regex_match(token,coef)||token=="+") 
					continue;
    			}
for(int l=0;l<N;l++)
std::cout<<product_coefs(k-1,l)<<"\t";
std::cout<<"\n";
		}
		else if(check=='b')
		{
			end =s.std::string::find("=",0);
std::cout<<"Back rxn\n";
		}
//to get a,n Ea for forward and backward rxns
		start =s.std::string::find_first_of("1234567890-+",end);
		end = s.std::string::find_first_of(" \t\n",start);
		token = s.substr(start, end - start);
		if(check=='b') Ar(k-1)=std::atof(token.c_str());
		else Af(k-1)=std::atof(token.c_str());

		start =s.std::string::find_first_of("1234567890-+",end);
		end = s.std::string::find_first_of(" \t\n",start);
		token = s.substr(start, end - start);
		if(check=='b') nr(k-1)=std::atof(token.c_str());
		else nf(k-1)=std::atof(token.c_str());

		start =s.std::string::find_first_of("1234567890-+",end);
		end = s.std::string::find_first_of(" \t\n",start);
		token = s.substr(start, end - start);
		if(check=='b') Ear(k-1)=std::atof(token.c_str());
		else Eaf(k-1)=std::atof(token.c_str());
std::cout<<Af(k-1)<<"\t"<<nf(k-1)<<"\t"<<Eaf(k-1)<<"\t"<<Ar(k-1)<<"\t"<<nr(k-1)<<"\t"<<Ear(k-1)<<"\n";


//to get pressure dependence coefs
		if((end)==std::string::npos)
		{
			std::getline(fin,s);
			start =s.std::string::find_first_of("1234567890-+",0);
			end = s.std::string::find_first_of(" \t\n",start);
			std::vector<double> Pcof;
			Pcof.push_back(k-1);
			while(start!=std::string::npos) 
			{
				token = s.substr(start, end - start);
				Pcof.push_back(std::atof(token.c_str()));
				start =s.std::string::find_first_of("1234567890-+",end);
				end = s.std::string::find_first_of(" \t\n",start);
			}
			std::getline(fin,s);
			start =s.std::string::find_first_of("1234567890-+",0);
			end = s.std::string::find_first_of(" \t\n",start);
			while(start!=std::string::npos) 
			{
				token = s.substr(start, end - start);
				Pcof.push_back(std::atof(token.c_str()));
				start =s.std::string::find_first_of("1234567890-+",end);
				end = s.std::string::find_first_of(" \t\n",start);
			}
std::cout<<"for:"<<k<<check<<"\n";
for(int l=0;l<10;++l) std::cout<<Pcof[l]<<"\t";
std::cout<<"\n";
Press_coefs.push_back(Pcof);
for(int l=0;l<10;++l) std::cout<<Press_coefs[0][l]<<"\t";
		}
		std::cout<<"\n";
	
	}
	fin.close();
	thirdbodies();
std::cout<<"done!\n";
}

void Properties::thirdbodies()
{
std::cout<<"in thirdbodies\n";
	//getting indices of rxns with third bodies..+ getting indices of Ms in Species
	std::vector<unsigned> rxni;
	std::regex tbody("M[0-9]+");
std::cout<<"M_index\n";
	for(unsigned i=0;i<N;++i)
	{
		if(std::regex_match(Species[i],tbody))
		{
			for(unsigned j=0;j<L;++j)
			{
				if(reactant_coefs(j,i)!=0)
					rxni.push_back(j);

			}
			M_index.push_back(rxni);
	//removing Ms from matrices
			Species.erase(Species.begin()+i);
			removeColumn(reactant_coefs,i);
			removeColumn(product_coefs,i);
			N--;i--;

		}
	}
std::cout<<"okay\n";

	//initialising Mxx and getting values from file
	Mxx.resize(M_index.size(),N); 
	Mxx=Eigen::MatrixXd::Zero(M_index.size(),N); 

	std::ifstream fin;
	fin.open(mechanismfile.c_str(),std::ios::in);
	std::regex Spc("[A-Z][A-Za-z0-9\\-]*");
	std::regex coef("[0-9\\.]*");
	std::regex Mspc("\\[[A-Z][A-Za-z0-9\\-]*\\]");
	std::regex Mline("[Let M[0-9]+[ A-Za-z0-9\\+=\\.\\[\\]]+");
	std::size_t start,end,old;
	std::string token,s;
	int mk,Spi;

	while(fin)
	{

		std::getline(fin,s);
		if(s.std::string::length()==0) continue;
		
		Spi=mk=0;
//!! Ms appear in the same order as in rxns
//Let M1 =  + 1.9 [CO] + 12 [H2O] + 2.5 [H2] + 3.8 [CO2] + 1.0 [OTHER].
		if(s.std::string::substr(0,5)=="Let M")
		{
			old=start=0;
			end = s.std::string::find_first_of(" \t\n",start);
			int lch=0;
			while(lch!=1)
    			{
				old=start;
				start = end + 1;
				end = s.std::string::find_first_of(" \t\n",start);
				if(end==std::string::npos)
				{
					end=s.length()-1;
					lch=1;
				}
				token = s.substr(start, end - start);
				if(std::regex_match(token,Mspc))
				{
					token=token.std::string::substr(1,end-start-2);
					Spi=get_Index(token,Species);
					if(std::regex_match(token,Spc)&&Spi!=-1)
					{
						std::string prev=s.std::string::substr(old,start-old-1);
						if(std::regex_match(prev,coef))
						{
							Mxx(mk,Spi)=std::atof(prev.c_str());
						}
						else
							Mxx(mk,Spi)=1;
					}
				}
				else 
					continue;
    			}
			mk++;
std::cout<<mk<<"\t";
for(unsigned i=0;i<N;++i)
std::cout<<Mxx(mk-1,i)<<"\t";
std::cout<<"\n";
		}
	}
std::cout<<"L="<<L<<"N="<<N<<"Spe size"<<Species.size()<<"\n";
}

