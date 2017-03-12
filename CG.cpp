//============================================================================
// Name        : CG.cpp
// Author      :
// Version     :
// Copyright   :
// Description : Defines the conjugate gradient method
//============================================================================
#include "waveguide.h"

//a function that returns a vector result from the multiplication of a nxn matrix and a nx1 vector
std::vector<double> Matrix_multiplier(std::vector<std::map<int, double> > A, std::vector<double> x)
{
	std::vector<double> A_x (A.size(),0.0);
	for (unsigned int i=0; i < A.size(); ++i)
    {
        for(std::map<int,double>::iterator it = A[i].begin(); it != A[i].end(); ++it)
        {
            A_x[i] += x[it->first]*it->second;
        }
    }
	return A_x;
}

//overloaded operator* to compute the dot product
double operator * (const std::vector<double> &vector1, const std::vector<double> &vector2)
{
    int i = 0;
    double result = 0;
    for (double first : vector1) {
        result += first * vector2[i];
        ++i;
    }
    return result;
}

std::vector<double> CG(std::vector<std::map<int, double> > A, std::vector<double> b)
{
	double oldres = 0.0;
    double newres = 0.0;
    double alpha;
    std::vector<double> r (A.size());
    std::vector<double> p (A.size());
    std::vector<double> x (A.size());
    std::vector<double> Ax (A.size());
    std::vector<double> Ap (A.size());
    
    //initialize
    Ax = Matrix_multiplier(A,x);
    for (unsigned int i=0; i < A.size(); ++i)
    {
        r[i] = b[i] - Ax[i];
    }
    p = r;
    oldres = r*r;
    
    //CG loop
    for (int k=0; k < 1e6; ++k)
    {
        //compute A*p
        Ap = Matrix_multiplier(A,p);
        
        //compute alpha
        alpha = oldres/(p*Ap);
        
        //compute x;
        for (unsigned int i=0; i < x.size(); ++i)
        {
            x[i] += alpha*p[i];
        }
        
        //compute r
        for (unsigned int i=0; i < r.size(); ++i)
        {
            r[i] -= alpha*Ap[i];
        }
        
        //compute new residual
        newres = r*r;
        if (sqrt(newres) < 1e-10)
        {
			break;
		}
		
        //compute new p
        for (unsigned int i=0; i < r.size(); ++i)
        {
            p[i] = r[i]+newres/oldres*p[i];
        }
        
        oldres = newres;
	}
    
    
    
    
	return x;
}

