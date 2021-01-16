#include "GenerateMu.h"

GenerateMu::GenerateMu(int a, std::string path_spectrum){
	interpolation_points=a;
	// в самый конец массива пишем значение плотности для данного элемента
	// (interpolation_points+1) to keep density value
	interpmu = new float[(interpolation_points+1)*2];
	interp_mu = new float*[(interpolation_points + 1)];
	for(int i = 0; i < (interpolation_points + 1); ++i){
		interp_mu[i] = interpmu + 2*i;}

	spectru_m = new float[interpolation_points * 2];
	s_pectrum = new float*[interpolation_points];
	for (int i = 0; i < interpolation_points; ++i) {
		s_pectrum[i] = spectru_m + 2 * i;
	}
	// парсер файла спектра источника
	FILE * spectrum_file;
	int row_counter = 0;
	char mystring[200];
	std::string read_out;
	total_hw = 0;
	fopen_s(&spectrum_file, path_spectrum.c_str(), "r");
	if (spectrum_file == NULL) perror("Error opening data file");
	else {
		while (!feof(spectrum_file)) {
			fgets(mystring, 200, spectrum_file);
			std::stringstream ss(mystring);
			float E_keV, p_hotons;
			ss >> E_keV >> p_hotons;
			if (p_hotons>0){
				s_pectrum[row_counter][0] = E_keV;
//				несущественная перенормировка от квадратных миллиметров к квадратным санитметрам
//				original spectrum is given in photons per mm^2 --> photons per cm^2
//				though visually it makes no difference in my algorithm
				s_pectrum[row_counter][1] = 1000*p_hotons;
				total_hw += p_hotons;
				row_counter++;}
		}
	}
	fclose(spectrum_file);
}

GenerateMu::~GenerateMu(){
// This destructor function for class Phantom deallocates the directory's list of Entry objects.
delete [] interpmu;
delete [] interp_mu;

delete[] spectru_m;
delete[] s_pectrum;
}

// парсер текстовых файлов базы данных по мю
void GenerateMu::mu_by_name(std::string f_ile, std::string d_ir){
	std::string t_xt= ".txt";
	std::string path_to_file=d_ir+f_ile+t_xt;
	std::cout<<path_to_file<<"\n";
	double rho=0;
	// Первый проход файла данных по мю используется только для информационных целей
	// в частности для подсчета строк данных row_counter
	// just to count the number of data rows in a NIST file
	int row_counter=0;
	FILE * data_file;
    char mystring [200];
	std::string read_out;
    fopen_s(&data_file, path_to_file.c_str(), "r");
		if (data_file == NULL) perror ("Error opening data file");
		else{while(!feof(data_file)){
				fgets (mystring ,200, data_file);
				
				if (strstr(mystring, "ro = ") != NULL){
					read_out = std::string(mystring);
					read_out.erase(0, 4);
					rho=atof(read_out.c_str());
					std::cout <<"element/compound density="<< rho << "g/cm^3"<<"\n";
				}
			// 27.12.2020 Here I detect lines with data by looking for the "E+" string
				if (strstr(mystring,"E+")!=NULL){
					row_counter++;
					std::stringstream ss(mystring);
					double E_keV, mu_1, mu2;
					// alltogether I read three (3) values from the line
					ss >> E_keV >> mu_1 >> mu2;
					//std::stringstream ss_new(mystring);
					// new data commented out on 07.01.2021
					// double E_Mev, Mu_1, Mu_2, Mu_3, Mu_4, Mu_5, Mu_6, Mu_7, Mu_8;
					// new data have more fields
					// ss_new >> E_Mev >> Mu_1 >> Mu_2 >> Mu_3 >> Mu_4 >> Mu_5 >> Mu_6 >> Mu_7;
					// std::cout << E_Mev <<"\t"<<Mu_1 << "\t" << Mu_2 << "\t" << Mu_3 << "\t" << Mu_4 << "\t" << Mu_5 << "\t" << Mu_6 << "\t" << Mu_7 << "\n";
				}
		}
		}
    fclose(data_file);
	
	// Пересчитать текстовый файл данных по мю из базы данных
	// для значений энергий присутствующих в спектре источника

	// Allocate dynamic array to store the original NIST data, i.e. energy and second column of the mu data
	int entries_count=row_counter;
	float* mudata = new float[entries_count*2];
	float** mu_data = new float*[entries_count];
	// initialize it
	for(int i = 0; i < entries_count; ++i){
		mu_data[i] = mudata + 2*i;}

	// теперь считываем данные в подготовленные массивы для дальнейшей интерполяции
	// согласно количеству записей в спектре источника
	// repeat the data redout again but now fill the data array with entries
	row_counter=0;
    fopen_s(&data_file,path_to_file.c_str(), "r");
		if (data_file == NULL) perror ("Error opening data file");
		else{while(!feof(data_file)){
				fgets (mystring ,200, data_file);
				if (strstr(mystring,"E+")!=NULL){
					std::stringstream ss(mystring);
					double E_keV, mu_1, mu2;
					ss >>E_keV>>mu_1>>mu2;
					mu_data[row_counter][0]= (float)E_keV;
					mu_data[row_counter][1]= (float)mu_1;
					row_counter++;}}}
    fclose(data_file);

	for(int i=0; i<interpolation_points;++i){ // цикл по ненулевым записям в спектре источника 
		interp_mu[i][0] = s_pectrum[i][0];
		double b_ottom= interp_mu[i][0];	double t_op = interp_mu[i][0];
		int f_lag=0; int index_bottom=0, index_top=0;

		for (int j=0;j<entries_count;++j){
			if( mu_data[j][0] <= interp_mu[i][0] ){b_ottom = mu_data[j][0]; index_bottom=j;}
			if( (mu_data[j][0]>= interp_mu[i][0])&&(f_lag!=1) ){t_op = mu_data[j][0];f_lag=1;index_top=j;break;}}
		
		// линейная интерполяция в логарифмических координатах
		// I'm using linear interpolation in log log space, i.e. connecting two data points by line

		interp_mu[i][1] = exp(log(mu_data[index_bottom][1]) + (log(mu_data[index_top][1]) - log(mu_data[index_bottom][1])) /
			(log(mu_data[index_top][0]) - log(mu_data[index_bottom][0]))*(log(interp_mu[i][0]) - log(mu_data[index_bottom][0])));

					 
		if(index_bottom==index_top){ interp_mu[i][1]=mu_data[index_bottom][1];};
		}
		interp_mu[interpolation_points][1] = (float)rho; // в конец массива пишем данные по плотности этого элемента
		delete[] mu_data;
		delete[] mudata;
}
