#include <iostream>
#include <stdlib.h>
#include "waveguide.h"

using namespace std;

double delta, epsilon;

int main(int argc, char** argv){


    int refLvl;
    if (argc < 3)
	{
		cout << "Input: (delta) (epsilon)" << endl;
		return -1;
	}
	else if(argc == 3)
    {
        delta = strtod(argv[1], NULL);
        epsilon = strtod(argv[2], NULL);
        refLvl = 0;
    }
	else if(argc == 4)
    {
        delta = strtod(argv[1], NULL);
        epsilon = strtod(argv[2], NULL);
        refLvl = atoi(argv[3]);
    }


    WaveGuide *wg = new WaveGuide(delta);

    wg->ReadInputFile();

    wg->Refine(refLvl);

    wg->CalculateKsqr();

    wg->WriteToFile("ksq.txt");

    wg->GetAhMhMatrices();
    
    wg->IPI();

	return 0;
}
