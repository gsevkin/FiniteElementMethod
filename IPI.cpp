//============================================================================
// Name        : IPI.cpp
// Author      :
// Version     :
// Copyright   :
// Description : Defines the Inverse Power Iteration function
//============================================================================

#include "waveguide.h"

extern double delta, epsilon;

//make a function to return the absolute value of input
double abs_d(double d)
{
	return (d>0.0) ? d:-d;
}

void WaveGuide::IPI()
{
    std::vector<double> f (_M.size());
	std::vector<double> u (_M.size(),1.0);
	double d = 1.0;
	double old_lambda;
	double lambda = 1.0;
	std::cout.precision(10);
	//(1)
	while(abs_d(d) > epsilon)
	{
		//(2)
		old_lambda = lambda;
		//(3)
		f = Matrix_multiplier(_M,u);
		//(4)
		u = CG(_A,f);
		//(5)
		double u_norm = sqrt(u*u);
		for (unsigned int i = 0; i < u.size(); ++i)
		{
			u[i] /= u_norm;
		}
		//(6)
		std::vector<double> Au = Matrix_multiplier(_A,u);
		std::vector<double> Mu = Matrix_multiplier(_M,u);
		lambda = (u*Au)/(u*Mu);
		std::cout << "Approximated lambda is: " << lambda << std::endl;
		//(7)
		d = (lambda - old_lambda)/old_lambda;
		
	}
	
	//print to file
	std::ofstream outFile;
    outFile.open("eigenmode.txt");

    for(int i = 0; i < num_vertex; ++i)
    {
        Coordinate c = _v.find(i)->second;
        outFile << c.x << ' ' <<c.y << ' ' << u[i] << std::endl;
    }

    outFile.close();
}
