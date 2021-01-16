#ifndef LOADPROJECTION_H
#define LOADPROJECTION_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <sstream>
#include <omp.h>
#include <math.h>

#include<algorithm>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

class Loadprojection {
	int y_source, y_detector;
	int w_detector, h_detector, projections_number, mammo;
	float start_angle, end_angle, cut_background, tukey_window_alpha;
	std::string read_projections_from;
	std::string a_ffix;
// конструктур и деструктор класса
public:
	int slice_size, t_arget, d_elta, times_n;
	Loadprojection( int, int, int, int,
					float, float,int, std::string, int, int,int,float,float,
					int, std::string, int);
	~Loadprojection();
	void load_images();
};
#endif
