#ifndef PHANTOM_H
#define PHANTOM_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <sstream>
#include <omp.h>

#include<algorithm>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

class Phantom {
	int w_idth, d_epth, h_eight, interpolation_points, elements_count, chamber_index;
	// Параметры размытия и зашумления проекции
	int gaussblur_kernel, gaussnoise_sigma, gaussnoise_mean;
	float a_multiplier;
	std::string write_projections_to;
	public:
		//  Немногочисленные переменные описывающие размеры фантома и его материал		
		uchar *ph_antom;

		Phantom(int, int, int, int);
		~Phantom();

		void create_phantom(uchar);
		void grow_random_defects(int, int, int, int);
		void save_phantom();
		void o_bject(int,int,int,int,int,int,int,float, float, float,int);
		void set_blur_and_noise(int,float,int,int);
		// Переменные используемые пр прорисовке проекции	
		int x_source, y_source, z_source;	
		int x_detector, y_detector, z_detector ;	
		void create_projection(int, int, int, float, int, int, float*[], float*[], int,int, float, int, int, float, float, int, std::string);
					
		};
#endif
