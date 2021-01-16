#ifndef GENERATEMU_H
#define GENERATEMU_H

#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <new>

class GenerateMu {
	public:
	int interpolation_points;
	// данные по мю персчитанные для спектра источника
	float* interpmu;
	float** interp_mu;
	// спектр источника
	float* spectru_m;
	float** s_pectrum;

	float total_hw;

	GenerateMu(int, std::string);
	~GenerateMu();
				
	void mu_by_name(std::string, std::string);
		};
		
#endif
