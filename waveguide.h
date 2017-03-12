/*
 * waveguide.h
 *
 *  Created on: Jun 9, 2015
 *      Author: uk74utuv
 */

#ifndef WAVEGUIDE_H_
#define WAVEGUIDE_H_

#include <string>
#include <math.h>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "Colsamm_Source/Colsamm.h"

using namespace ::_COLSAMM_;


struct Coordinate
{
	double x;
	double y;

	bool operator <( const Coordinate &rhs ) const
    {
       return ( x < rhs.x);
    }

};

struct Face
{
	int v0;
	int v1;
	int v2;
};

class WaveGuide {
public:

	WaveGuide(double delta);
	virtual ~WaveGuide();

	void ReadInputFile();

	void CalculateKsqr();

	void WriteToFile(std::string fileName);

    void GetAhMhMatrices();

    void Refine(int refLvl);

    void IPI();

private:

	std::map<int, Coordinate> _v;
	std::vector<Face> _f;
    std::map<int, double> _ksqr;
    std::vector<std::map<int, double>> _A;
    std::vector<std::map<int, double>> _M;

    int num_vertex = 1039;
	int num_face = 1976;

    void RefineOnce();
    void setNFace(int);
    void setNVertex(int);
};

std::vector<double> Matrix_multiplier(std::vector<std::map<int, double> > A, std::vector<double> x);
double operator * (const std::vector<double> &vector1, const std::vector<double> &vector2);
std::vector<double> CG(std::vector<std::map<int, double> > A, std::vector<double> b);

#endif /* WAVEGUIDE_H_ */
