class Mechanism {
	read_ThermoData();
	//	Input: Thermodynamic data file
	//	Output: Species names, Mass, heat of formation, entropy of formation, heat capacity
	read_Mechanism();
	//	Input: Mechanism file
	//	Output variables: Species names, Reactions
	//	Reactions will involves reactant species, product species and coefficients, rate constant : A, n, Ea
	read_Transport();
	//	Input: Transport file
	//	Output: Collision diameter, viscosity, thermal conductivity, collision energy
};

class ZeroD {
	
};