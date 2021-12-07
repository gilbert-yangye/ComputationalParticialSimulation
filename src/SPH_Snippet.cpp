#include "SPH_2D.h"
#include "file_writer.h"
#include <omp.h>
#include <cmath>
#include <fstream>

SPH_main domain;

int main(void)
{
	// Domain configuration
	double length = 20.0, height = 10.0, water_depth = 2.0;
	double src_length = 3.0, src_height = 3.0, dx = 0.2;
	int domain_case = 1;

	double start = omp_get_wtime();  // start timing

	domain.set_values(length, height, dx);	// set simulation parameters

	domain.initialise_grid();	// initialise simulation grid

	// places initial points - include the boundary points and the specifics of where the fluid is in the domain
	// for case 2 and case 3, can change the steepness in the function implementation
	domain.initialise_domain(water_depth, src_length, src_height, domain_case);				


	// Time-stepping configuration
	double final_t = 15.0;
	double dt = 0.1 * domain.h / domain.c_0;
	double c_cfl = 0.1;

	// for post-processing in python
	ofstream fname;
	fname.open("dt.txt");
	fname << dt << endl;
	fname.close();

	domain.time_stepping(dt, final_t, c_cfl, domain_case);

	double end = omp_get_wtime(); // end timing
	cout << "Running time = " << end - start << endl;

	return 0;
}
