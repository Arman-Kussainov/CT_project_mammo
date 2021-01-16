#include "Phantom.h"

using namespace cv;
using namespace std;

const float pi = 3.14159f;
const float rad_grad = 0.01745f;

Phantom::Phantom(int a, int b, int c, int chamberid) {
	w_idth = a; d_epth = b; h_eight = c; chamber_index = chamberid;
	ph_antom = new (std::nothrow) uchar[w_idth*d_epth*h_eight];
	if (ph_antom == 0) { std::cout << "\n Error assigning memory \n"; }
}

Phantom::~Phantom(){
// деструктор класса 
delete [] ph_antom;
}

// заполняем фантом его материалом
void Phantom::create_phantom(uchar phantom_material) {
	std::fill(ph_antom, ph_antom + w_idth*d_epth*h_eight, phantom_material);}
// "выращивание" дефектов случайных по форме, размеру и расположению 
// тип _1 соответсвует "непустым" дефектам "defects", тип _2 соответсвует "пустотам"-"voids"
void Phantom::grow_random_defects(int material_1,int quantity_1,
					   int material_2,int quantity_2){
	int defect_material[2], defect_quantity[2];
	defect_material[0]=material_1; defect_material[1]=material_2;
	defect_quantity[0]=quantity_1; defect_quantity[1]=quantity_2;
							   			
	srand(static_cast<unsigned char>(time(NULL)));
// цикл по "дефектам" и "пустотам"
	for (int d_efects=0;d_efects<=1;++d_efects){
// цикл по заявленному количеству "дефектов"
		for(int ie=1;ie<=defect_quantity[d_efects];ie++){
// случайным образом выбирается тип дефектов, эллипсы или прямоугольники
			int defect_type=rand()%2;
// максимальные и минимальные размеры дефектов по отношению к фантому
			int max_size=w_idth/4;
			int min_size=w_idth/8;
// пересчет размеров самих дефектов или пустот						
			int widt_h= min_size + rand()%max_size;
			int dept_h= min_size + rand()%max_size;
			int heigh_t=min_size + rand()%max_size;
// положение объекта внутри фантома
			float x= (float)(rand()%w_idth);
			float y= (float)(rand()%d_epth);
			float z= (float)(rand()%h_eight);
// Подобранные случайным образом углы поворота объектов
			float alpha= float(rand())/RAND_MAX*2.0f*pi;
			float beta = float(rand())/RAND_MAX*2.0f*pi;
			float gamma= float(rand())/RAND_MAX*2.0f*pi;
// заранее опсчитанные 3D матрицы поворота точки объекта
			float Rx[3][3]={1,0,0,
							 0,cos(alpha),-sin(alpha),
							 0,sin(alpha),cos(alpha)};
			float Ry[3][3]={cos(beta),0,sin(beta),
							 0,1,0,
							-sin(beta),0,cos(beta)};
			float Rz[3][3]={cos(gamma),-sin(gamma),0,
							sin(gamma),cos(gamma),0,
							0,0,1};
// заполняем и поворачиваем на заданные углы все точки текущего объекта (дефекта или пустоты)
			for(int k=0;k<widt_h;++k){
				for(int l=0;l<dept_h;++l){
					for(int m=0;m<heigh_t;++m){

						int half_object_width = widt_h / 2;
						int half_object_depth = dept_h / 2;
						int half_object_height = heigh_t / 2;

						float ko = static_cast<float>(k - half_object_width);
						float lo = static_cast<float>(l - half_object_depth);
						float mo = static_cast<float>(m - half_object_height);


						// поворот вокруг оси Х
						float k1 = ko*Rx[0][0] + lo*Rx[0][1] + mo*Rx[0][2];
						float l1 = ko*Rx[1][0] + lo*Rx[1][1] + mo*Rx[1][2];
						float m1 = ko*Rx[2][0] + lo*Rx[2][1] + mo*Rx[2][2];
						// поворот вокруг оси У
						float k2 = k1*Ry[0][0] + l1*Ry[0][1] + m1*Ry[0][2];
						float l2 = k1*Ry[1][0] + l1*Ry[1][1] + m1*Ry[1][2];
						float m2 = k1*Ry[2][0] + l1*Ry[2][1] + m1*Ry[2][2];
						// поворот вокруг оси Z + смещение + и переход от значений координат к значениям индексов 3D матрицы
						// за счет static_cast<int>
						int X= static_cast<int>(k2*Rz[0][0]+l2*Rz[0][1]+m2*Rz[0][2]+x);
						int Y= static_cast<int>(k2*Rz[1][0]+l2*Rz[1][1]+m2*Rz[1][2]+y);
						int Z= static_cast<int>(k2*Rz[2][0]+l2*Rz[2][1]+m2*Rz[2][2]+z);
						// заполнение объекта его материалом если после всех преобразований и поворотов он оказался внутри 
						// фантома
						// Эти три цикла пришлось ввести из-за ошибки дискретизации и округления при повороте фантома
						// при этом повернутый объект становился нозреватым и незаполненным полностью материалом фантома
						for (int Xi = X - 1; Xi <= X + 1;Xi++){
							for (int Yi = Y - 1; Yi <= Y + 1; Yi++){
								for (int Zi = Z - 1; Zi <= Z + 1; Zi++) {

						if((Xi>=0)&&(Yi>= 0)&&(Zi>=0)&&(Xi<=w_idth-1)&&(Yi<=d_epth-1)&&(Zi<=h_eight-1)&&
						// не забываем, что пустоты-ваакум съедают все
						// однако если фантом не из "вакуума" и мы начинаем выращивать строго "непустые" дефекты
						// такая преодсторожность избыточна...
						(ph_antom[Xi + Yi*w_idth + Zi*w_idth*d_epth]!=0)){
							if(defect_type==0){// эллипсоиды
								if( pow(ko, 2.0f) / (float)(half_object_width*half_object_width) +
									pow(lo, 2.0f) / (float)(half_object_depth*half_object_depth) +
									pow(mo, 2.0f) / (float)(half_object_height*half_object_height) <= 1.0f)
									{ph_antom[Xi + Yi*w_idth + Zi*w_idth*d_epth]=defect_material[d_efects];}}
							else{ //параллелепипеды
								ph_antom[Xi + Yi*w_idth + Zi*w_idth*d_epth]=defect_material[d_efects];}}
								
						}}}

			}}}
		}
	}
}

// все тоже самое только для детерминированных объектов
// исходные английские комментарии сохранены
void Phantom::o_bject(int object_shape,
					  int object_x,int object_y,int object_z,
					  int object_width,int object_depth,int object_height,
	float object_alpha, float object_beta, float object_gamma,
					  int object_id){			  

	float Rx[3][3]={ 1,0,0,
					 0,cos(object_alpha),-sin(object_alpha),
					 0,sin(object_alpha),cos(object_alpha)};
	float Ry[3][3]={cos(object_beta),0,sin(object_beta),
					 0,1,0,
					-sin(object_beta),0,cos(object_beta)};
	float Rz[3][3]={cos(object_gamma),-sin(object_gamma),0,
					sin(object_gamma),cos(object_gamma),0,
					0,0,1};
	for(int k=0;k<object_width;++k){
		for(int l=0;l<object_depth;++l){
			for(int m=0;m<object_height;++m){

				int half_object_width =  object_width  / 2;
				int half_object_depth =  object_depth  / 2;
				int half_object_height = object_height / 2;

				float ko = static_cast<float>(k- half_object_width);
				float lo = static_cast<float>(l - half_object_depth);
				float mo = static_cast<float>(m - half_object_height);

				// rotation around X axis
				float k1=ko*Rx[0][0]+lo*Rx[0][1]+mo*Rx[0][2];
				float l1=ko*Rx[1][0]+lo*Rx[1][1]+mo*Rx[1][2];
				float m1=ko*Rx[2][0]+lo*Rx[2][1]+mo*Rx[2][2];
				// rotation around Y axis
				float k2=k1*Ry[0][0]+l1*Ry[0][1]+m1*Ry[0][2];
				float l2=k1*Ry[1][0]+l1*Ry[1][1]+m1*Ry[1][2];
				float m2=k1*Ry[2][0]+l1*Ry[2][1]+m1*Ry[2][2];
				// rotation around z axis + translation + transform to index notation
				int X= static_cast<int>(k2*Rz[0][0]+l2*Rz[0][1]+m2*Rz[0][2]+float(object_x));
				int Y= static_cast<int>(k2*Rz[1][0]+l2*Rz[1][1]+m2*Rz[1][2]+float(object_y));
				int Z = static_cast<int>(k2*Rz[2][0] + l2*Rz[2][1] + m2*Rz[2][2] + float(object_z));
				// growing the defect inside the phantom's matrix

				// Эти три цикла пришлось ввести из-за ошибки дискретизации и округления при повороте фантома
				// без них, повернутый объект становился ноздреватым и незаполненным полностью материалом фантома
				for (int Xi = X - 1; Xi <= X + 1; Xi++) {
					for (int Yi = Y - 1; Yi <= Y + 1; Yi++) {
						for (int Zi = Z - 1; Zi <= Z + 1; Zi++) {
							if ((Xi >= 0) && (Yi >= 0) && (Zi >= 0) && (Xi <= w_idth - 1) && (Yi <= d_epth - 1) && (Zi <= h_eight - 1) &&
								// next line makes sure that voids destroys everything
								(ph_antom[Xi + Yi*w_idth + Zi*w_idth*d_epth]) != 0) {
								if (object_shape == 0) {// ellipsoids
									if (pow(ko, 2.0f) / (float)(half_object_width*half_object_width) +
										pow(lo, 2.0f) / (float)(half_object_depth*half_object_depth) +
										pow(mo, 2.0f) / (float)(half_object_height*half_object_height) <= 1.0f)
									{
										ph_antom[Xi + Yi*w_idth + Zi*w_idth*d_epth] = object_id;
									}
								}
								else { //cuboids
									ph_antom[Xi + Yi*w_idth + Zi*w_idth*d_epth] = object_id;
								}
							}

				}}}


	}}}					  
}

// сохранение фантома в виде текстового файла на диск
void Phantom::save_phantom(){
// первые три элемента матрицы содержат информацию о размерах фантома!
// это было необходимо для визуализации в Матлаб
	double dtime1 = omp_get_wtime();
	FILE* fout;
	fopen_s(&fout, "phantom.txt", "wb");
	for (int k = 0; k < h_eight; ++k) {
		for (int j = 0; j < d_epth; ++j) {
			for (int i = 0; i < w_idth; ++i) {
				int c_ounter = i + j*w_idth + k*w_idth*d_epth;
				if (c_ounter == 0) { fprintf(fout, "%u\t", w_idth); }
		   else if (c_ounter == 1) { fprintf(fout, "%u\t", d_epth); }
		   else if (c_ounter == 2) { fprintf(fout, "%u\t", h_eight); }
				else { fprintf(fout, "%u\t", ph_antom[c_ounter]); }
			}fprintf(fout, "\n");}}
	fclose(fout);
	dtime1 = omp_get_wtime() - dtime1;
	std::printf("time to save phantom to text file in seconds = %f\n\n", dtime1);
}

// функция задающая парметры размытых и зашумленных проекций
void Phantom::set_blur_and_noise(int gaussblurkernel, float a_value, int noisemean, int gaussnoisesigma) {
	gaussblur_kernel = gaussblurkernel;
	gaussnoise_sigma = gaussnoisesigma;
	gaussnoise_mean=noisemean;
	a_multiplier =a_value;}

// основная процедура трассировки объекта и записи проекции при заданном угле "theta_phantom"
void Phantom::create_projection(int y_source, int z_source, int y_detector,
								float theta_phantom,
								int w_detector, int h_detector,
								float** interpolated_mu, float** s_pectrum , int interpolation_points, int elements_count,
								float total_photons, int image_compression, int voxels_per_cm, float t_before, float t_after, int mammo,
								std::string projections_path){
	write_projections_to = projections_path;
	float a_ngle = theta_phantom;
	
	t_after >= 1 ? t_after = 1 : t_after = t_after;
	cout << "y_source" << "\t" << "z_source" << "\t" << "y_detector" << "\t" << "t_before" << "\t" << "t_after" << "\n";
	cout << y_source << "\t\t" << z_source << "\t\t" << y_detector << "\t\t" << t_before << "\t" << t_after << "\n";
	
	double dtime = omp_get_wtime();
	// максимальное число параметров описывающие переход от вращающейся системы координат привязанной к 
	// источнику/детектору к неподвижной фантома вынесено за пределы вложенных циклов для того чтобы увеличить скорость
	// вычислений. Например вычисление синусов и косинусов углов поворота записывающей системы
	float c_os=cos(theta_phantom*rad_grad);
	float s_in=sin(theta_phantom*rad_grad);

	float y_dist = float(y_detector - y_source); 
	float y_dist_cos = y_dist*c_os;
	float y_dist_sin = y_dist*s_in;

	float y_source_sin = float(y_source)*s_in - float(w_idth) / 2.0f;
	float y_source_cos = float(y_source)*c_os + float(d_epth) / 2.0f;
	float z_source_add = float(z_source) + float(h_eight) / 2.0f;

	float x_width_half = float(w_idth) / 2.0f;
	float y_d_epth_half = float(y_source) + float(d_epth) / 2.0f;
	float z_h_eight_half = float(z_source) + float(h_eight) / 2.0f;

	float voxel_size = 1.0f / float(voxels_per_cm); // cm

// Если источник находится не на оси проходящей через центр объекта и параллельной его (объекта) оси ОУ
// вдоль которой ведется съемка, а например выше ее, нам нужно смещать детектор вниз, чтобы
// держать в поле зрения тень объекта
	float z_detector = float(z_source*y_detector) / float(abs(y_source));
	// включение выключение режима маммографии
	if (mammo == 1) { z_detector = 0; }

	float total_flux = 0;
// сумарное количество фотонов у источника
// нужно для однородной нормировки всех проекций, для различных углов
	for (int k = 0; k < interpolation_points; k++) {
		total_flux += s_pectrum[k][1];
	}
// ИНТЕРЕСНЫЙ ФАКТ! ПРИ УВЕЛИЧЕНИИ РАЗМЕРОВ ФАНТОМА И ДЕТЕКТОРА, СВЕТОСИЛА УСТАНОВКИ РАСТЕТ!
// А ИМЕННО, КОЛИЧЕСТВО ТРАССИРУЕМЫХ ЛУЧЕЙ РАВНО КОЛИЧЕСТВУ ПИКСЕЛЕЙ В ДЕТЕКТОРЕ!
// ТО ЕСТЬ, ЧЕМ БОЛЬШЕ ДЕТЕКТОР, ТЕМ БОЛЬШЕ ЛУЧЕЙ ГЕНЕРИРУЕТ ИСТОЧНИК

	Mat img_float = cv::Mat(h_detector, w_detector, CV_32FC1, Scalar::all(0)); // исходная проекция/изображение
	Mat img;  // исходная проекция/изображение

	Mat g_auss; // изображение размытое гауссовской функцией
	Mat n_oise = cv::Mat(h_detector, w_detector, CV_16UC1); // изображение с белым шумом

	Mat s_um; // для sqrt(1+a*img_float) данных
	Mat noise_float; // белый шум type float
	
	// другие вспомогательные Mat объекты
	Mat add_noise = cv::Mat(h_detector, w_detector, CV_32FC1);
	Mat img_noise = cv::Mat(h_detector, w_detector, CV_32FC1);
	Mat img_noise_uchar = cv::Mat(h_detector, w_detector, CV_16UC1);



// Инструкции для многоядерных процессоров использовать 
// максимальное количество ядер и соответственно потоков им доступное
// Explicitly disable dynamic teams	
omp_set_dynamic(0);     
// Get the number of processors in this system
int c_ores = omp_get_num_procs();
// Now set the number of threads
omp_set_num_threads(c_ores);
// иногда необходимо протестировать один только поток thread
//omp_set_num_threads(1);

#pragma omp parallel 
	{
// количество вокселей какого-либо типа на пути луча
		float * voxel_encounters;
		voxel_encounters = new (std::nothrow) float[elements_count + 1]();
		if (voxel_encounters == 0) { std::cout << "\n Error assigning memory \n"; }
// этот спектр мы будем модифицировать для каждого из лучей, при прохождении чере вещество
		float * modified_photons;
		modified_photons = new (std::nothrow) float[interpolation_points]();
		if (modified_photons == 0) { std::cout << "\n Error assigning memory \n"; }
// заранее пересчитывет значение экспоненциального ослабления для каждого из элементов на шаге в один воксель
		float *e_xp;
		e_xp = new (std::nothrow) float[interpolation_points*(elements_count + 1)]();
		if (e_xp == 0) { std::cout << "\n Error assigning memory \n"; }

		for (int i = 0; i < interpolation_points; i++) {
			for (int j = 1; j <= elements_count; j++) { // j=0 соответствует вакууму и не используется.
				// элемент interpolated_mu[j][interpolation_points] содержит значение плотности данного вещества rho
				// dt равно одному вокселю во всех направлениях, см. ниже
				e_xp[j + i*elements_count] = exp(-interpolated_mu[j][i] * interpolated_mu[j][interpolation_points] * voxel_size);
			}
		}

// Это блок предназначен для корректного описания фантома Шеппа-Логана через относительные плотности веществ
// [1, 0.2, 0, 0, 0.3, 0.4, 0.3, 0.3, 0.4, 0]
// ***	
//		float rho[11] = { 0.001f, 1.0f, 0.2f, 0.0f, 0.0f, 0.3f, 0.4f, 0.3f, 0.3f, 0.4f, 0.0f };
//		for (int i = 0; i < interpolation_points; i++) {
//			for (int j = 1; j <= elements_count; j++) { // j=0 соответствует вакууму и не используется.
//														// элемент interpolated_mu[j][interpolation_points] содержит значение плотности данного вещества rho
//														// dt равно одному вокселю во всех направлениях, см. ниже
//				e_xp[j + i*elements_count] =
//					exp(-interpolated_mu[2][i]*interpolated_mu[2][interpolation_points] * voxel_size * rho[j-1]);
//			}
//		}
// ***


#pragma omp for
		for (int j = 0; j<img_float.rows; j++) {
			
			float z_dist = float(j) - float(h_detector) / 2.0f - z_detector - float(z_source);

			for (int i = 0; i<img_float.cols; i++) {
			float i_ntensity = 0;

// 12.05.2018 перенесено и объеденино с другим циклом ниже поо тексту
// создаем копию спектра у источника и будем ее модифицировать для каждого из лучей
//			for (int k = 0; k < interpolation_points; k++) {
//				modified_photons[k] = s_pectrum[k][1];}

			float x_dist = float(i) - float(w_detector) / 2.0f;
			float xcos_ysin = x_dist*c_os - y_dist_sin;
			float xsin_ycos = x_dist*s_in + y_dist_cos;

// Если фантом находится в камере заполненной веществом, надо учитывать это при трассировке луча
// сначала ослабляем каждый элемент спектра на весь путь в данном веществе
// потом, при трассировке отменяем ослабление для каждого индивидуального вокселя заполненного веществом фантома
			float dist_in_cm = pow(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist, 0.5f)* voxel_size;

				for (int energy_bins = 0; energy_bins < interpolation_points; energy_bins++) {
// создаем копию спектра у источника и будем ее модифицировать для каждого из лучей
					modified_photons[energy_bins] = s_pectrum[energy_bins][1];
					if (chamber_index > 0) {
						modified_photons[energy_bins] *=
							float(exp(-interpolated_mu[chamber_index][energy_bins] * interpolated_mu[chamber_index][interpolation_points]
							* dist_in_cm));}						
				}

			
			float dt = pow((y_dist*y_dist + x_dist*x_dist + z_dist*z_dist), -0.5f);
			float inv_r2 = 1/(y_dist*y_dist + x_dist*x_dist + z_dist*z_dist);
			// процедура трасировки индивидуального луча для всех элементов спектра			
			int Xt, Yt, Zt;
			for(float t= t_before;t<=t_after;t+=dt){ // tracing only through the region occupied by a phantom
				if (mammo == 0) {
					Xt = static_cast<int>(xcos_ysin*t - y_source_sin); //
					Yt = static_cast<int>(xsin_ycos*t + y_source_cos);
					Zt = static_cast<int>(z_dist*t + z_source_add);
				}
				else {
			// в режиме маммографии фантом находится в неопсредственной близости от детектора
					Xt = static_cast<int>( x_dist*t + x_width_half);
					Yt = static_cast<int>( y_dist*t + y_d_epth_half);
					Zt = static_cast<int>( z_dist*t + z_h_eight_half);
				}

			// производим модификацию спектра луча в зависимости от материала фантома в данной точке
			// при этом надо обратить экспоненциальное ослабление луча в этой точке коэфицентом мю вещества камеры 
			// проведенное РАНЕЕ для всего луча, в том числе внутри фантома, см. цикл по спектру источника выше
				if ((Xt >= 0) && (Yt >= 0) && (Zt >= 0) && (Xt < w_idth) && (Yt < d_epth) && (Zt < h_eight)) {
					int inde_x = ph_antom[Xt + Yt*w_idth + Zt*w_idth*d_epth];
					
					for (int energy_bins = 0; energy_bins < interpolation_points; energy_bins++) {
						
						if (chamber_index > 0) { // испытательная камера заполнена каким-нибудь материалом
							float reverse_air_voxel = e_xp[chamber_index + energy_bins*elements_count];
							modified_photons[energy_bins] /= reverse_air_voxel; // обращаем ослабление материалом камеры
							if (inde_x != 0) { // ослабляем поток фотонов материалом вокселя
								modified_photons[energy_bins] *= e_xp[inde_x + energy_bins*elements_count]; }}
						// если материал испытательной камеры вакуум, просто ослабляем поток фотонов материалом вокселя
						else if (inde_x != 0) {modified_photons[energy_bins] *= e_xp[inde_x + energy_bins*elements_count];}
					}
				}

			
			}

			// собираем весь спектр (кол-во фотонов на участок спектра) в одну интенсивность
			for (int energy_bins = 0; energy_bins < interpolation_points; energy_bins++) {
				i_ntensity += modified_photons[energy_bins];}
			// можно моделировать монохроматических спектр
			//i_ntensity +=modified_photons[25];}
			
			// ослабляем луч обратно пропорциональнео квадрату расстояния			
			i_ntensity *= inv_r2;
			img_float.at<float>(h_detector - 1 - j, w_detector - 1 - i) = i_ntensity;
	}}
	
		
	}

// массивы задаваемыми динамическим способом необходимо удлять из памяти
// при работе с многими параллельными потоками заново их декларируя 
// в параллельной среде
#pragma omp parallel 
	{
		float * modified_photons;
		modified_photons = new (std::nothrow) float[interpolation_points]();
		if (modified_photons == 0) { std::cout << "\n Error assigning memory \n"; }
		delete[]modified_photons;

		float *e_xp;
		e_xp = new (std::nothrow) float[interpolation_points*(elements_count + 1)]();
		if (e_xp == 0) { std::cout << "\n Error assigning memory \n"; }
		delete[]e_xp;

	}

	// нормировка всех проекций на интенсивность центрального луча
	// ослабленного на квадрат расстояния и не встретившего на пути никаких препятствий
	// Тем самым мы инвертируем изображение!
	total_flux /=(y_dist*y_dist);
	img_float = abs(img_float - total_flux);
	cv::normalize(img_float, img_float, 0, 65535, CV_MINMAX);
	img_float.convertTo(img, CV_16UC1);

	vector<int> compression_params;
	compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
	compression_params.push_back(image_compression);

	string pa_th = write_projections_to;

	string projectio_n="projection_";
	string extensio_n=".png";
	string extension_gauss= "_gauss.png";
	string extension_rand = "_rand.png";
	string extension_noise = "_noise.png";
	
	string counte_r = SSTR(a_ngle);
		
	string filenam_e = pa_th+projectio_n + counte_r + extensio_n;
	string filename_gauss = pa_th + projectio_n + counte_r + extension_gauss;
	string filename_noise = pa_th + projectio_n + counte_r + extension_noise;

	std::cout << filenam_e<<"\n";
	// запись размытых и зашумленных проекций не проводится в данной версии программы
	// std::cout << filename_gauss << "\n";
	// std::cout << filename_noise << "\n";

	// Размытие и зашумление проекций согласно заданию Проекта  № 1
	// 1. Gaussian blur
	cv::GaussianBlur(img, g_auss, Size(gaussblur_kernel, gaussblur_kernel), 0, 0);
	// 2. Adding white noise
	// another uchar to float conversion, i.e. working with smaller 0-255 range of values now
	img.convertTo(img_float, CV_32FC1);
	float a = a_multiplier;
	cv::sqrt(1+a*img_float, s_um); // s_um is float
	
	cv::randn(n_oise, gaussnoise_mean, gaussnoise_sigma); // unsigned char like in original projection
	n_oise.convertTo(noise_float, CV_32FC1); // conversion without scaling
	cv::multiply(noise_float,s_um, add_noise); // multiplication
	
	img_noise = img_float + add_noise; //  + шум
	img_noise.convertTo(img_noise_uchar, CV_16UC1); // clipping above 65535 is done automaticaly 
	// запись размытых и зашумленных проекций не проводится в данной версии программы
 	cv::imwrite(filenam_e.c_str(), img, compression_params); 
	// cv::imwrite(filename_gauss.c_str(), g_auss, compression_params);
	// cv::imwrite(filename_noise.c_str(), img_noise_uchar, compression_params);

	dtime = omp_get_wtime() - dtime;
	std::printf("elapsed build time in seconds = %f\n\n", dtime);
}
