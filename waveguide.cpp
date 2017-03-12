//============================================================================
// Name        : waveguide.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Defines all the class member function
//============================================================================
#include "waveguide.h"

extern double delta;

WaveGuide::WaveGuide(double delta)
{
    
    _A.resize(num_vertex);
    _M.resize(num_vertex);
}

WaveGuide::~WaveGuide() {
	// TODO Auto-generated destructor stub
}

void WaveGuide::Refine(int refLvl)
{
    for(int i=0; i<refLvl; ++i)
    {
        RefineOnce();
    }
}
Coordinate FindMidPoint(Coordinate x1, Coordinate x2)
{
    Coordinate c;
    c.x = (x1.x+x2.x)/2;
    c.y = (x1.y+x2.y)/2;
    return c;
}

void WaveGuide::RefineOnce()
{
    int vert_index=1;
    std::map<int, Coordinate> nv;
	std::vector<Face> nf;

    std::map<Coordinate, int> v_ind;

    for(int i=0; i< num_face; ++i)
    {
        // new vertex
        Coordinate c1 = _v.find(_f[i].v0)->second;
        Coordinate c2 = _v.find(_f[i].v1)->second;
        Coordinate c3 = _v.find(_f[i].v2)->second;

        Coordinate m1 = FindMidPoint(c1, c2);
        Coordinate m2 = FindMidPoint(c1, c3);
        Coordinate m3 = FindMidPoint(c2, c3);

        // add vertex index
        if ( v_ind.find(c1) == v_ind.end() )
        {
            v_ind.insert(std::pair<Coordinate, int>(c1, vert_index));
            nv.insert(std::pair<int, Coordinate>(vert_index, c1));
            ++vert_index;
        }
        if ( v_ind.find(c2) == v_ind.end() )
        {
            v_ind.insert(std::pair<Coordinate, int>(c2, vert_index));
            nv.insert(std::pair<int, Coordinate>(vert_index, c2));
            ++vert_index;
        }
        if ( v_ind.find(c3) == v_ind.end() )
        {
            v_ind.insert(std::pair<Coordinate, int>(c3, vert_index));
            nv.insert(std::pair<int, Coordinate>(vert_index, c3));
            ++vert_index;
        }
        if ( v_ind.find(m1) == v_ind.end() )
        {
            v_ind.insert(std::pair<Coordinate, int>(m1, vert_index));
            nv.insert(std::pair<int, Coordinate>(vert_index, m1));
            ++vert_index;
        }
        if ( v_ind.find(m2) == v_ind.end() )
        {
            v_ind.insert(std::pair<Coordinate, int>(m2, vert_index));
            nv.insert(std::pair<int, Coordinate>(vert_index, m2));
            ++vert_index;
        }
        if ( v_ind.find(m3) == v_ind.end() )
        {
            v_ind.insert(std::pair<Coordinate, int>(m3, vert_index));
            nv.insert(std::pair<int, Coordinate>(vert_index, m3));
            ++vert_index;
        }

        Face tmp_f;
        // new face 1
        tmp_f.v0 = v_ind.find(c1)->second;
        tmp_f.v1 = v_ind.find(m1)->second;
        tmp_f.v2 = v_ind.find(m2)->second;
        nf.push_back(tmp_f);

        // new face 2
        tmp_f.v0 = v_ind.find(c2)->second;
        tmp_f.v1 = v_ind.find(m1)->second;
        tmp_f.v2 = v_ind.find(m3)->second;
        nf.push_back(tmp_f);

        // new face 3
        tmp_f.v0 = v_ind.find(c3)->second;
        tmp_f.v1 = v_ind.find(m2)->second;
        tmp_f.v2 = v_ind.find(m3)->second;
        nf.push_back(tmp_f);

        // new face 4
        tmp_f.v0 = v_ind.find(m1)->second;
        tmp_f.v1 = v_ind.find(m2)->second;
        tmp_f.v2 = v_ind.find(m3)->second;
        nf.push_back(tmp_f);
    }

    // reset local variables
    _f.clear();
    _f = nf;
    setNFace(4*num_face);

    _v.clear();
    _v = nv;
    setNVertex(vert_index);

    _A.resize(num_vertex);
    _M.resize(num_vertex);
}

void WaveGuide::setNFace(int nf)
{
    num_face = nf;
}

void WaveGuide::setNVertex(int nv)
{
    num_vertex = nv;
}

double KsqrFunction(double x, double y)
{
    double k2 = (100 + delta)*exp(-50.0*(x*x + y*y)) -100;
    return k2;
}

void WaveGuide::GetAhMhMatrices()
{
    ELEMENTS::Triangle my_element;
    std::vector< std::vector< double > > my_local_matrix;
    std::vector<double> corners(6, 0.0);

    for(int i = 0; i < num_face; ++i)
    {
        Coordinate c1 = _v.find(_f[i].v0)->second;
        Coordinate c2 = _v.find(_f[i].v1)->second;
        Coordinate c3 = _v.find(_f[i].v2)->second;

        corners[0] = c1.x;
        corners[1] = c1.y;
        corners[2] = c2.x;
        corners[3] = c2.y;
        corners[4] = c3.x;
        corners[5] = c3.y;

        // pass the corners to the finite element
        my_element(corners);

        my_local_matrix = my_element.integrate(grad(v_()) * grad(w_()) - func<double>(KsqrFunction) * v_() * w_());

        _A[_f[i].v0][_f[i].v0] += my_local_matrix[0][0];
        _A[_f[i].v0][_f[i].v1] += my_local_matrix[0][1];
        _A[_f[i].v0][_f[i].v2] += my_local_matrix[0][2];

        _A[_f[i].v1][_f[i].v0] += my_local_matrix[1][0];
        _A[_f[i].v1][_f[i].v1] += my_local_matrix[1][1];
        _A[_f[i].v1][_f[i].v2] += my_local_matrix[1][2];

        _A[_f[i].v2][_f[i].v0] += my_local_matrix[2][0];
        _A[_f[i].v2][_f[i].v1] += my_local_matrix[2][1];
        _A[_f[i].v2][_f[i].v2] += my_local_matrix[2][2];

        my_local_matrix = my_element.integrate(v_() * w_());

        _M[_f[i].v0][_f[i].v0] += my_local_matrix[0][0];
        _M[_f[i].v0][_f[i].v1] += my_local_matrix[0][1];
        _M[_f[i].v0][_f[i].v2] += my_local_matrix[0][2];

        _M[_f[i].v1][_f[i].v0] += my_local_matrix[1][0];
        _M[_f[i].v1][_f[i].v1] += my_local_matrix[1][1];
        _M[_f[i].v1][_f[i].v2] += my_local_matrix[1][2];

        _M[_f[i].v2][_f[i].v0] += my_local_matrix[2][0];
        _M[_f[i].v2][_f[i].v1] += my_local_matrix[2][1];
        _M[_f[i].v2][_f[i].v2] += my_local_matrix[2][2];

    }

    // print A matrix
    std::ofstream outFileA, outFileM;
    outFileA.open("A.txt");
    outFileM.open("M.txt");
    for(int i = 0; i < num_vertex; ++i)
    {
        for(std::map<int,double>::iterator itA = _A[i].begin(); itA != _A[i].end(); ++itA)
        {
            outFileA << i << ' ' << itA->first << ' ' << itA->second << std::endl;
        }

        for(std::map<int,double>::iterator itM = _M[i].begin(); itM != _M[i].end(); ++itM)
        {
            outFileM << i << ' ' << itM->first << ' ' << itM->second << std::endl;
        }

    }
    outFileA.close();
    outFileM.close();
}

void WaveGuide::WriteToFile( std::string fileName )
{
    std::ofstream outFile;
    outFile.open(fileName);

    for(int i = 0; i < num_vertex; ++i)
    {
        Coordinate c = _v.find(i)->second;
        outFile << c.x << "\t" <<c.y << "\t" << _ksqr.find(i)->second << std::endl;
    }

    outFile.close();
}


void WaveGuide::CalculateKsqr()
{
    for(int i = 0; i < num_vertex; ++i)
    {
        Coordinate c = _v.find(i)->second;
        double k2 = (100 + delta)*exp(-50.0*(c.x*c.x + c.y*c.y)) -100;

        _ksqr.insert(std::pair<int, double>(i, k2));
    }
}

void WaveGuide::ReadInputFile()
{

	std::ifstream infile("./inputs/unit_circle.txt");
	std::string line;
	int line_counter = 0;

	while (std::getline(infile, line))
	{
		std::istringstream iss(line);

	    if(line_counter < num_vertex)
	    {
		    int ver_num;
		    Coordinate v;
		    if (!(iss >> ver_num >> v.x >> v.y)) { continue; }
		    ++line_counter;
		    _v.insert ( std::pair<int,Coordinate>(ver_num,v) );
	    }

	    else if(line_counter< num_vertex + num_face )
	    {
            Face f;
   		    if (!(iss >> f.v0 >> f.v1 >> f.v2)) { continue; }
		    ++line_counter;
            _f.push_back(f);
	    }

        else
        {
            break;
        }
	}

}
