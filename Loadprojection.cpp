#include "Loadprojection.h"

using namespace cv;
using namespace std;

const float pi = 3.14159f;
const float rad_grad = 0.01745f;

Loadprojection::Loadprojection(int ysource, int ydetector, int wdetector, int hdetector,
	float startangle, float endangle, int F_OV, string slice_plane, int target_cross_section, int slice_half_thickness, int zero_padding_factor,
	float cutbackground, float tukeywindowalpha,
	int projectionsnumber, string projections_path, int m_ammo) {
	y_source = ysource;
	y_detector = ydetector;
	w_detector = wdetector; h_detector = hdetector;
	projections_number = projectionsnumber;
	slice_size = F_OV;
	a_ffix = slice_plane;
	t_arget = target_cross_section;
	d_elta = slice_half_thickness;
	times_n = zero_padding_factor;
	cut_background = cutbackground;
	tukey_window_alpha = tukeywindowalpha;
	read_projections_from = projections_path;
	start_angle = startangle;
	end_angle = endangle;
	mammo = m_ammo;
}

Loadprojection::~Loadprojection() {}
// Функция перестановки левой и правой части матриц

Mat swapLR(Mat complex)
{
	Mat tmp, q0, q1;
	// Мы работаем с матрицами имеющими четное количество столбцов!
	int cx = complex.cols / 2;
	int cy = complex.rows;
	// Копируем левые и правые части наших матриц в промежуточные переменные
	q0 = complex(Rect(0, 0, cx, cy));
	q1 = complex(Rect(cx, 0, cx, cy));
	// Меняем скопированные части местами
	q0.copyTo(tmp);
	q1.copyTo(q0);
	tmp.copyTo(q1);
	return complex;
}
// Фурье преобразование над массивом одномерных данных
Mat compute_rowDFT(Mat image) {
	Mat planes[] = { Mat_< float>(image), Mat::zeros(image.size(), CV_32F) };
	// для комплексных чисел, нужно два набора данных Re и Im
	Mat complex;
	merge(planes, 2, complex);
	// Фурье преобразование 
	dft(complex, complex, DFT_ROWS | DFT_COMPLEX_OUTPUT);
	return complex;
}

// функция "дозабивания" нулями
Mat pad_zeros(Mat image, int zeros_pad) {
	Mat padded;
	int m = 0, n = 0, d_ff = zeros_pad - image.cols;
	if (d_ff % 2 == 0) {
		m = d_ff / 2; n = d_ff / 2;
	}
	else { m = d_ff / 2; n = m + 1; }
	// create output image of optimal size
	//copyMakeBorder( src, dst, top, bottom, left, right, borderType, value );
	copyMakeBorder(image, padded, 0, 0, m, n, BORDER_CONSTANT, Scalar::all(0));
	return padded;
}

Mat crop_zeros(Mat image, int fit_size) {
	int m = 0, n = 0, d_ff = image.cols - fit_size;
	if (d_ff % 2 == 0) {
		m = d_ff / 2; n = d_ff / 2;
	}
	else { m = d_ff / 2; n = m + 1; }
	Mat cropped = image(Rect(m, 0, image.rows, image.rows)).clone();
	return cropped;
}

// От обратного преобразования Фурье нам нужна только действительная часть planes[0]
Mat inverse_rowDFT(Mat complex)
{ // CHECK OUT THE SCALLING!!!!
	Mat work;
	// normalization for inverse Fourier transform
	cv::divide(complex, complex.cols, complex);
	dft(complex, work, DFT_ROWS | DFT_INVERSE | DFT_COMPLEX_OUTPUT);
	Mat planes[] = { Mat::zeros(complex.size(), CV_32F), Mat::zeros(complex.size(), CV_32F) };
	split(work, planes);                // planes[0] = Re(DFT(I)), planes[1] = Im(DFT(I))
	//magnitude(planes[0], planes[1], work);    // === sqrt(Re(DFT(I))^2 + Im(DFT(I))^2)
	//return work;
	return planes[0];
}

// Произведение двух комплексных матриц
Mat mult_complx(Mat complex1, Mat complex2)
{
	Mat planes1[] = { Mat::zeros(complex1.size(), CV_32F), Mat::zeros(complex1.size(), CV_32F) };
	split(complex1, planes1);                // planes[0] = Re(DFT(I)), planes[1] = Im(DFT(I))

	Mat planes2[] = { Mat::zeros(complex2.size(), CV_32F), Mat::zeros(complex2.size(), CV_32F) };
	split(complex2, planes2);                // planes[0] = Re(DFT(I)), planes[1] = Im(DFT(I))

	Mat Re_Im[] = { Mat::zeros(complex1.size(), CV_32F), Mat::zeros(complex1.size(), CV_32F) };
	Re_Im[0] = planes1[0].mul(planes2[0]) - planes1[1].mul(planes2[1]);
	Re_Im[1] = planes1[0].mul(planes2[1]) + planes1[1].mul(planes2[0]);
	Mat p_rod;
	merge(Re_Im, 2, p_rod);
	return p_rod;
}

void Loadprojection::load_images() {
	double dtime = omp_get_wtime();

	vector<int> compression_params;
	compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
	compression_params.push_back(0);

	float view_angle;
	float half_slice = float(slice_size) / 2.0f;

	Mat c_onvolved = cv::Mat(h_detector, w_detector, CV_32FC1, Scalar::all(0));

	Mat reconstructed_image = cv::Mat(slice_size, slice_size, CV_32FC1, Scalar::all(0));
	Mat reconstructed_scaled;

	Mat IMG_float = cv::Mat(h_detector, w_detector, CV_32FC1, Scalar::all(0));

	Mat gna_filter = cv::Mat(h_detector, w_detector*times_n, CV_32FC1, Scalar::all(0));
	Mat hann_filter = cv::Mat(h_detector, w_detector*times_n, CV_32FC1, Scalar::all(0));
	Mat Tukey_filter = cv::Mat(h_detector, w_detector, CV_32FC1, Scalar::all(0));

	Mat IMG = cv::Mat(h_detector, w_detector, CV_16UC1, Scalar::all(0));

	Mat slice_normalized = cv::Mat(slice_size, slice_size, CV_16UC1, Scalar::all(0));
	Mat slice_all = cv::Mat(slice_size, slice_size, CV_32FC1, Scalar::all(0));

	float *g_na;
	float *h_ann;
	float *T_ukey;

	// Если мы увеличиваем размер проекции нулями, тоже самое нужно сделать с фильтрами
	// times_n это целочисленный фактор на который мы увеличиваем фильтры в пространстве частот
	int filter_width = w_detector*times_n;
	g_na = new (std::nothrow) float[filter_width]();
	h_ann = new (std::nothrow) float[filter_width]();
	T_ukey = new (std::nothrow) float[w_detector]();

	// 1D Tukey window Окно Таки в пространстве исходной проекции
	// see https://en.wikipedia.org/wiki/Window_function
	float a_lpha = tukey_window_alpha;
	if (a_lpha > 0) {
		int min_dim = w_detector;
		if (h_detector < w_detector) { min_dim = h_detector; }

		int range_1 = static_cast<int>(a_lpha*static_cast<float>(min_dim - 1) / 2.0f);
		int range_2 = static_cast<int>(static_cast<float>(min_dim - 1)*(1.0f - a_lpha / 2.0f));

		for (int ni = 0; ni < min_dim; ni++) {
			if ((0 <= ni) && (ni < range_1)) {
				T_ukey[ni] = 0.5f*(1.0f + cos(pi*(2.0f*static_cast<float>(ni) / (a_lpha*static_cast<float>(min_dim - 1)) - 1.0f)));
			}

			if ((range_1 <= ni) && (ni <= range_2)) {
				T_ukey[ni] = 1.0;
			}

			if ((range_2 < ni) && (ni <= min_dim - 1)) {
				T_ukey[ni] = 0.5f*(1.0f + cos(pi*(2.0f*static_cast<float>(ni) / (a_lpha*static_cast<float>(min_dim - 1)) - 2.0f / a_lpha + 1.0f)));
			}

		}

		int center_ck = static_cast<int>(float(h_detector) / 2.0f);
		int center_ci = static_cast<int>(float(w_detector) / 2.0f);

		for (int ck = 0; ck < h_detector; ck++) {
			for (int ci = 0; ci < w_detector; ci++) {
				int dist_to_cntr = static_cast<int>(
					pow((float)((ci - center_ci)*(ci - center_ci) + (ck - center_ck)*(ck - center_ck)), 0.5));

				int half_dist = static_cast<int>(float(min_dim) / 2.0f);
				if (dist_to_cntr <= half_dist) {
					Tukey_filter.at<float>(ck, ci) = T_ukey[half_dist - dist_to_cntr];
				}
			}
		}
	}
	else { Tukey_filter = Scalar::all(1); }

	// Считаем значения элементов фильтра g_na в действительном пространстве
	// а фильтр h_ann в пространстве частот Фурье
	int center_point = static_cast<int>(float(filter_width) / 2.0f);
	int A = 1;
	g_na[center_point] = 1 / float(8 * A*A);

	for (int ni = 1; ni <= filter_width / 2; ni++) {
		int index_left = center_point - ni;
		int index_right = center_point + ni;
		if (ni % 2 != 0) {
			if (index_left >= 0) { g_na[index_left] = -1 / (2 * float(ni*ni*A*A)*pi*pi); }
			if (index_right < filter_width) { g_na[index_right] = -1 / (2 * float(ni*ni*A*A)*pi*pi); }
		}
	}

	// окно Хамминга 
	for (int ni = 0; ni < filter_width; ni++) {
		h_ann[ni] = 0.54f - 0.46f*cos(2.0f*pi*float(ni) / float(filter_width - 1));
	}

	// заполняем 1D фильтром 2D матрицу
	for (int ck = 0; ck < h_detector; ck++) {
		for (int ci = 0; ci < w_detector*times_n; ci++) {
			gna_filter.at<float>(h_detector - 1 - ck, w_detector*times_n - 1 - ci)
				= g_na[w_detector*times_n - 1 - ci];
			hann_filter.at<float>(h_detector - 1 - ck, w_detector*times_n - 1 - ci)
				= h_ann[w_detector*times_n - 1 - ci];
		}
	}

	// переводим наш фильтр в пространство частот
	Mat gna_fft = compute_rowDFT(gna_filter);
	// делаем двухслойным для перемножения на двухслойный комплексный gna_fft
	// чтобы перемножить на gna_fft поэлементно, так как hann_filter сугубо действительный
	hann_filter = swapLR(hann_filter);
	Mat two_hann[] = { Mat::zeros(hann_filter.size(), CV_32F), Mat::zeros(hann_filter.size(), CV_32F) };
	two_hann[0] = hann_filter;
	two_hann[1] = hann_filter;
	Mat ha_nn;
	merge(two_hann, 2, ha_nn);
	gna_fft = gna_fft.mul(ha_nn);

	// фактор нормирующий изначальную проекцию, см. отчет
	Mat Dna = cv::Mat(h_detector, w_detector, CV_32FC1, Scalar::all(0));
	float dist_to_screen = float(-y_source + y_detector);
	for (int ck = 0; ck < h_detector; ck++) {
		float z_dist = float(ck) - float(h_detector) / 2.0f;
		for (int ci = 0; ci < w_detector; ci++) {
			float x_dist = float(ci) - float(w_detector) / 2.0f;
			Dna.at<float>(ck, ci) =
				dist_to_screen / pow(dist_to_screen*dist_to_screen + x_dist*x_dist + z_dist*z_dist, 0.5f);
		}
	}

	string extensio_n = ".png";
	//	string extensio_n = ".tif";

	// Вычисления шага поворота пары источник-детектор
	float degree_step;
	if (projections_number == 1) {
		// любое значение отличное от нуля, чтобы цикл выполнялся один раз, для одной проекции
		degree_step = 1.0;
	}
	else {
		degree_step = (end_angle - start_angle) / (float(projections_number) - 1.0f);
	}

	for (float d_egree = start_angle; d_egree <= end_angle; d_egree += degree_step) {
		if (mammo == 0) {
			view_angle = (d_egree - start_angle) * pi / 180;
		}
		else { view_angle = d_egree * pi / 180; }

		// вычисляем тригонометрический функции углов поворота за пределами циклов
		float sin_view_angle = sin(view_angle);
		float cos_view_angle = cos(view_angle);

		// фактор нормирующий модифицированную проекцию, см. отчет
		Mat U2 = cv::Mat(slice_size, slice_size, CV_32FC1, Scalar::all(0));
		// calculate U2 weighting
#pragma omp parallel for
		for (int x_t = 0; x_t < slice_size; x_t++) {
			for (int y_t = 0; y_t < slice_size; y_t++) {
				U2.at<float>(x_t, y_t) =
					(float)pow(float(-y_source) /
					(float(-y_source) + (float(x_t) - half_slice)*sin_view_angle + (float(y_t) - half_slice)*cos_view_angle), 2.0);
			}
		}

		// Восстанавливаемое сечение или FOV
		Mat slice = cv::Mat(slice_size, slice_size, CV_32FC1, Scalar::all(0));

		float dist_to_screen = float(-y_source + y_detector);

		// Считывание и подготовка проекций для обработки
		string counte_r = SSTR(d_egree);

		string projectio_n = "projection_";
		//		string projectio_n = "RetrvdImg_";
		string p_ath = read_projections_from;

		string filenam_e = p_ath + projectio_n + counte_r + extensio_n;
		cout << "reading data from: " << filenam_e << "\n";

		IMG = imread(filenam_e, CV_LOAD_IMAGE_UNCHANGED);

		// Проверка наличия файла проекции
		if (!IMG.data) {
			cout << "Could not open or find the image" << std::endl;
		}
		IMG.convertTo(IMG_float, CV_32FC1);

		// Переход от интенсивностей к распределению коэфицентов мю, через логарифм
		IMG_float = IMG_float * 35.0f + 9e6f;
		log(IMG_float, IMG_float);
		cv::normalize(IMG_float, IMG_float, 0, 65535, CV_MINMAX);

		// Внимание zerro-padding должно применяться только к изображению,
		// чья интенсивность спадает к нулю по краям изображения
		IMG_float = IMG_float.mul(Tukey_filter);

		// Применение первого модифицирующего множителя, см.отчет		
		IMG_float = IMG_float.mul(Dna);
		Mat IMG_padded;
		if (times_n > 1) { IMG_padded = pad_zeros(IMG_float, w_detector*times_n); }
		else { IMG_padded = IMG_float; }
		/*------------------*/

// Переход в пространство Фурье для применения фильтра
		Mat complexI = compute_rowDFT(IMG_padded);
		/*------------------*/
		c_onvolved = inverse_rowDFT(mult_complx(complexI, gna_fft));
		//c_onvolved = inverse_rowDFT(complexI);
		c_onvolved = swapLR(c_onvolved);
		// возвращаем изображение увеличенное нулями в исходное состояние
		if (times_n > 1) { c_onvolved = crop_zeros(c_onvolved, w_detector); }

		// Поворачиваем восстанавливаемый объект в обратную сторону!  см. отчет
		float c_os = cos(-view_angle);
		float s_in = sin(-view_angle);
		float y_source_float = float(-y_source);

		// Демонстрация режима маммографии (не совсем корректна с математической точки зрения, потому что нет фильтров)
		if ((mammo == 1) && (a_ffix == "YZ")) {

			float source_to_detector_center = static_cast<float>(-y_source + y_detector);
			// в случае маммографии, старые переменные приобретают новый смысл и меняются с каждым шагом по углу
			// возвышения источника. Также, для маммографии, нам нужны положительные углы s_in -> (-s_in)
			// и знак их снова меняется с минуса на плюс, см. выше 19.01.2021
			float sourcetophantom = source_to_detector_center*c_os - static_cast<float>(y_detector);
			float sourceelevation = source_to_detector_center*(-s_in);
			float source_to_detector = source_to_detector_center*c_os; //cosine is an even function, so the sign of angle is irrelevant

			// распараллеливание процесса обратного проецирования
			for (int Xt = t_arget - d_elta; Xt <= t_arget + d_elta; Xt++) {
				float Xp = static_cast<float>(Xt) - half_slice;
				//if (abs((int)Xp) < y_detector) {
#pragma omp parallel for
					for (int Yt = 0; Yt < slice_size; Yt++) {
						// Yt and Zt are the coordinates on the reconstructed slice
						float Yp = static_cast<float>(Yt) - half_slice;
						// cut the area of the slice beyond the available space limited by a detector moved immediately to the phantom
						//if (abs((int)Yp) < y_detector) {
							int ci = static_cast<int>(Xp*source_to_detector / (sourcetophantom + Yp) + static_cast<float>(w_detector) / 2.0);
							// ci and ck are the coordinates of the detector 
							for (int Zt = 0; Zt < slice_size; Zt++) {
								float Zp = static_cast<float>(Zt) - half_slice;
								//if (abs((int)Zp) < y_detector) {
									int ck = static_cast<int>(sourceelevation - (sourceelevation - Zp)*source_to_detector / (sourcetophantom + Yp)
										+ static_cast<float>(h_detector) / 2.0);

									if ((ck >= 0) && (ck < h_detector) && (ci >= 0) && (ci < w_detector)) {
										slice.at<float>(Zt, Yt) += c_onvolved.at<float>(h_detector - 1 - ck, w_detector - 1 - ci);
									}

								//}
							}
						//}
					}
				//}
			}
		}

		// January 18, 2021. Trying to reproduce the old solution for the "YZ" mammo plane with the other planes

		if ((mammo == 1) && (a_ffix == "XY")) {

			float source_to_detector_center = static_cast<float>(-y_source + y_detector);
			// в случае маммографии, старые переменные приобретают новый смысл и меняются с каждым шагом по углу
			// возвышения источника. Также, для маммографии, нам нужны положительные углы s_in -> (-s_in)
			float sourcetophantom = source_to_detector_center*c_os - static_cast<float>(y_detector);
			float sourceelevation = source_to_detector_center*(-s_in);
			float source_to_detector = source_to_detector_center*cos(view_angle);


			for (int Zt = t_arget - d_elta; Zt <= t_arget + d_elta; Zt++) {
				float Zp = static_cast<float>(Zt - half_slice);
				//if (abs((int)Zp) < y_detector) {
					// распараллеливание процесса обратного проецирования
#pragma omp parallel for
					for (int Xt = 0; Xt < slice_size; Xt++) {
						float Xp = static_cast<float>(Xt) - half_slice;
						//if (abs((int)Xp) < y_detector) {
							for (int Yt = 0; Yt < slice_size; Yt++) {
								float Yp = static_cast<float>(Yt) - half_slice;
								// cut the area of the slice beyond the available space limited by a detector moved immediately to the phantom
								//if (abs((int)Yp) < y_detector) {
									int ci = static_cast<int>(Xp*source_to_detector / (sourcetophantom + Yp)
										+ static_cast<float>(w_detector) / 2.0);
									int ck = static_cast<int>(sourceelevation - (sourceelevation - Zp)*source_to_detector / (sourcetophantom + Yp)
										+ static_cast<float>(h_detector) / 2.0);

									if ((ck >= 0) && (ck < h_detector) && (ci >= 0) && (ci < w_detector)) {
										slice.at<float>(Xt, Yt) += c_onvolved.at<float>(h_detector - 1 - ck, w_detector - 1 - ci);
									}
								//}
							}
						//}
					}
				//}
			}
			// Для сечения параллельного плоскости XY, модифицирующий фактор можно вынести за циклы
			//slice = slice.mul(U2);
		}

		if ((mammo == 1) && (a_ffix == "ZX")) {

			float source_to_detector_center = static_cast<float>(-y_source + y_detector);
			float sourcetophantom = source_to_detector_center*c_os - static_cast<float>(y_detector);
			float sourceelevation = source_to_detector_center*(-s_in);
			float source_to_detector = source_to_detector_center*cos(view_angle);

			for (int Yt = t_arget - d_elta; Yt <= t_arget + d_elta; Yt++) {
				float Yp = static_cast<float>(Yt) - half_slice;
				// cut the area of the slice beyond the available space limited by a detector moved immediately to the phantom
				//if (abs((int)Yp) < 120) {

#pragma omp parallel for
					for (int Xt = 0; Xt < slice_size; Xt++) {
						float Xp = static_cast<float>(Xt) - half_slice;
						//if (abs((int)Xp) < 120) {
							int ci = static_cast<int>(Xp*source_to_detector / (sourcetophantom + Yp)
								+ static_cast<float>(w_detector) / 2.0);

							for (int Zt = 0; Zt < slice_size; Zt++) {
								float Zp = static_cast<float>(Zt) - half_slice;
								//if (abs((int)Zp) < 120) {
									int ck = static_cast<int>(sourceelevation - (sourceelevation - Zp)*source_to_detector / (sourcetophantom + Yp)
										+ static_cast<float>(h_detector) / 2.0);

									if ((ck >= 0) && (ck < h_detector) && (ci >= 0) && (ci < w_detector)) {
										slice.at<float>(Zt, Xt) += (c_onvolved.at<float>(h_detector - 1 - ck, w_detector - 1 - ci));
									}

								//}
							}
						//}
					}
				//}
			}
		}


		// для трех различных типов проекций используем свой блок комманд
		// ограничивающих выборку координат X,Y,Z с их переводом в координаты
		// проецируемых пикселей ci и ck рентгеновских проекций
		if ((mammo == 0) && (a_ffix == "XY")) {

			for (int Zt = t_arget - d_elta; Zt <= t_arget + d_elta; Zt++) {
				float Zp = static_cast<float>(Zt - half_slice);
				// распараллеливание процесса обратного проецирования
#pragma omp parallel for
				for (int Xt = 0; Xt < slice_size; Xt++) {
					float Xt_cos = static_cast<float>(Xt - half_slice)*c_os;
					float Xt_sin = static_cast<float>(Xt - half_slice)*s_in;
					for (int Yt = 0; Yt < slice_size; Yt++) {
						float Y_t = static_cast<float>(Yt - half_slice);

						float Xp = Xt_cos - Y_t*s_in;
						float Yp = Xt_sin + Y_t*c_os;

						float r_atio = dist_to_screen / (y_source_float + Yp);

						int ci = static_cast<int>(Xp*r_atio + static_cast<float>(w_detector) / 2.0);
						int ck = static_cast<int>(Zp*r_atio + static_cast<float>(h_detector) / 2.0);

						if ((ck >= 0) && (ck < h_detector) && (ci >= 0) && (ci < w_detector)) {
							slice.at<float>(Xt, Yt) += c_onvolved.at<float>(h_detector - 1 - ck, w_detector - 1 - ci);
						}


					}
				}
			}
			// Для сечения параллельного плоскости XY, модифицирующий фактор можно вынести за циклы
			slice = slice.mul(U2);
		}
		// все то же самое но для других перепендикулярных сечений
		if ((mammo == 0) && (a_ffix == "YZ")) {

			for (int Xt = t_arget - d_elta; Xt <= t_arget + d_elta; Xt++) {
				float Xt_cos = static_cast<float>(Xt - half_slice)*c_os;
				float Xt_sin = static_cast<float>(Xt - half_slice)*s_in;
#pragma omp parallel for
				for (int Yt = 0; Yt < slice_size; Yt++) {
					float Y_t = static_cast<float>(Yt - half_slice);
					// from rotated to fixed frame or vice versa
					float Xp = Xt_cos - Y_t*s_in;
					float Yp = Xt_sin + Y_t*c_os;

					float r_atio = dist_to_screen / (y_source_float + Yp);

					int ci = static_cast<int>(Xp*r_atio + static_cast<float>(w_detector) / 2.0);

					for (int Zt = 0; Zt < slice_size; Zt++) {

						float Zp = static_cast<float>(Zt - half_slice);
						int ck = static_cast<int>(Zp*r_atio
							+ static_cast<float>(h_detector) / 2.0);

						//check backward readout once more time
						if ((ck >= 0) && (ck < h_detector) && (ci >= 0) && (ci < w_detector)) {
							slice.at<float>(Zt, Yt) += (c_onvolved.at<float>(h_detector - 1 - ck, w_detector - 1 - ci)*U2.at<float>(Xt, Yt));
						}

					}
				}
			}

		}

		if ((mammo == 0) && (a_ffix == "ZX")) {

			for (int Yt = t_arget - d_elta; Yt <= t_arget + d_elta; Yt++) {
				float Yt_sin = static_cast<float>(Yt - half_slice)*s_in;
				float Yt_cos = static_cast<float>(Yt - half_slice)*c_os;
#pragma omp parallel for
				for (int Xt = 0; Xt < slice_size; Xt++) {
					float X_t = static_cast<float>(Xt - half_slice);

					float Xp = X_t*c_os - Yt_sin;
					float Yp = X_t*s_in + Yt_cos;

					float r_atio = dist_to_screen / (y_source_float + Yp);

					int ci = static_cast<int>(Xp*r_atio + static_cast<float>(w_detector) / 2.0);

					for (int Zt = 0; Zt < slice_size; Zt++) {

						float Zp = static_cast<float>(Zt - half_slice);
						int ck = static_cast<int>(Zp*r_atio + static_cast<float>(h_detector) / 2.0);


						if ((ck >= 0) && (ck < h_detector) && (ci >= 0) && (ci < w_detector)) {
							slice.at<float>(Zt, Xt) += (c_onvolved.at<float>(h_detector - 1 - ck, w_detector - 1 - ci)*U2.at<float>(Xt, Yt));
						}

					}
				}
			}
		}

		// Складываем информацию от каждой проекции в общий файл
		slice_all = slice_all + slice;

		cv::normalize(slice_all, slice_normalized, 0, 65535, CV_MINMAX, CV_16UC1);
		cv::namedWindow(a_ffix, WINDOW_NORMAL);
		cv::imshow(a_ffix, slice_normalized);
		cv::waitKey(1);
	}

	// создание и применение маски М для удаления фона изоборажения
	double min, max;
	cv::minMaxLoc(slice_all, &min, &max);
	min = min*(1.0f - cut_background);
	min = 400.0;
	Mat M_ask;
	inRange(slice_all, min, max, M_ask);
	M_ask.convertTo(M_ask, CV_32FC1, 1.f / 255);
	//	slice_all=slice_all.mul(M_ask);
	//	slice_all = slice_all / 3000.0 * 0.2;
	//	slice_all.setTo(1.0, slice_all > .7);

	cv::normalize(slice_all, slice_normalized, 0, 65535, CV_MINMAX, CV_16UC1);

	//	обработка исходного изображения и запись на диск
	string n_ame = "_2018_Dna_U2_";
	//string pat_h = "C:\\Users\\Martha\\Desktop\\slices\\";
	string pat_h = "C:\\Users\\Martha\\Desktop\\mammo_data\\";
	//string pat_h = "";
	string c_ounter = SSTR(t_arget);
	string time_stamp = SSTR(round(omp_get_wtime()));
	string rotation_s = SSTR(projections_number);
	string slic_e = SSTR(slice_size);
	string filenam_e = pat_h + c_ounter +"_"+rotation_s + "_" + a_ffix + "_"+slic_e + extensio_n;
	//string filenam_e = pat_h + c_ounter + "_" + a_ffix + n_ame + rotation_s + "_" + slic_e + "_" + time_stamp + extensio_n;
	cv::imwrite(filenam_e, slice_normalized, compression_params);

	dtime = omp_get_wtime() - dtime;
	std::printf("elapsed build time in seconds = %f\n\n", dtime);
	cv::imshow(a_ffix, slice_normalized);
	cv::waitKey(0);

	delete[]g_na;
	delete[]h_ann;
	delete[]T_ukey;
}
