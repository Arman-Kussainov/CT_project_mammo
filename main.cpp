#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <sstream>
#include <omp.h>

#include "Phantom.h"
#include "GenerateMu.h"
#include "Loadprojection.h"

const float pi = 3.14159f;
using namespace std;

int main() {
	double main_time = omp_get_wtime();

	int tet_variable = 44;

	// parse input file for the phantom info
	FILE * data_file;
	FILE * spectrum_file;
	string d_ir, spectrum_dir, projections_dir, slice_plane;

	char mystring[200] = {};
	string read_out;

	// последовательное считывание информации из управляющего файла input_file.txt
	// reading info about width, depth, height of the phantom and other info
	int width = 0, depth = 0, height = 0,
		phantom_id = 0, chamber_id = 0, defects_quantity = 0, defects_id = 0, voids_quantity = 0,
		source_to_phantom = 0, source_elevation = 0, phantom_to_detector = 0,
		F_OV = 0, target_cross_section = 0, slice_half_thickness = 0, zero_padding_factor = 1,
		detector_width = 0, detector_height = 0,
		projections_number = 0,
		image_compression = 0, voxels_per_cm = 0,
		gauss_blur_kernel = 0, gauss_noise_mean = 0, gauss_noise_sigma = 0,
		save_phantom_yes = 0, mammography_yes = 0, populate_projections = 0;
	float start_angle = 0.0f, end_angle = 0.0f, a_multiplier = 0.0f, cut_background=0.0f, tukey_window_alpha = 0.0f;

	// переменная используется для информации о количестве ненулевых значений количества фотонов
	// во всем спектре в предоставленном файле
	int interpolation_points = 0;

	// значение материала пустоты, далее по тексту, можно сделать ненулевым, для исследовательских целей
	// default int voids_id=0 can be overloaded by nonzero value from an input_file.txt	
	int voids_id = 0;
	fopen_s(&data_file, "input_file.txt", "r");
	if (data_file == NULL) perror("Error opening data file");
	else {
		while (!feof(data_file)) {
			fgets(mystring, 200, data_file);

			if (strstr(mystring, "path.mu.data") != NULL) {
				d_ir = string(mystring); d_ir.erase(0, 13);
				size_t found = d_ir.find_first_of("\n");
				d_ir.erase(int(found), 100);
				cout << "mu data files are at: " << d_ir << "\n";
			}

			if (strstr(mystring, "path.projection.data") != NULL) {
				projections_dir = string(mystring); projections_dir.erase(0, 21);
				size_t found = projections_dir.find_first_of("\n");
				projections_dir.erase(int(found), 100);
				cout << "projections files are at: " << projections_dir << "\n";
			}

			if (strstr(mystring, "path.spectrum.data") != NULL) {
				spectrum_dir = string(mystring); spectrum_dir.erase(0, 19);
				size_t f_ound = spectrum_dir.find_first_of("\n");
				spectrum_dir.erase(int(f_ound), 100);
				cout << "spectrum data file is: " << spectrum_dir << "\n";

				char storage_string[200];
				fopen_s(&spectrum_file, spectrum_dir.c_str(), "r");
				if (spectrum_file == NULL) perror("Error opening data file");
				else {
					while (!feof(spectrum_file)) {
						fgets(storage_string, 200, spectrum_file);
						double E_keV, p_hotons;
						stringstream ss(storage_string);
						ss >> E_keV >> p_hotons;
						if (p_hotons > 0) { interpolation_points++; }
					}
				}
				fclose(spectrum_file);
				cout << interpolation_points << " non zero data points are found in spectrum\n";
			}

			if (strstr(mystring, "source.to.phantom") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						source_to_phantom = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nsource.to.phantom = " << source_to_phantom;
			}

			if (strstr(mystring, "source.elevation") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						source_elevation = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nsource.elevation = " << source_elevation;
			}

			if (strstr(mystring, "phantom.to.detector") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						phantom_to_detector = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nphantom.to.detector = " << phantom_to_detector << "\n";
			}

			if (strstr(mystring, "field.of.view") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						F_OV = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nfield.of.view = " << F_OV << "\n";
			}

			if (strstr(mystring, "slice.plane") != NULL) {
				slice_plane = string(mystring); slice_plane.erase(0, 12);
				size_t found = slice_plane.find_first_of("\n");
				slice_plane.erase(int(found), 100);
				cout << "slice plane is: " << slice_plane << "\n";
			}

			if (strstr(mystring, "target.cross.section") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						target_cross_section = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "target.cross.section = " << target_cross_section << "\n";
			}

			if (strstr(mystring, "slice.half.thickness") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						slice_half_thickness = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "slice.half.thickness = " << slice_half_thickness << "\n";
			}

			if (strstr(mystring, "zero.padding.factor") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						zero_padding_factor = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "zero.padding.factor = " << zero_padding_factor;
			}

			if (strstr(mystring, "cut.background") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atof(read_out.c_str()) != 0) {
						cut_background = (float)atof(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ncut.background = " << cut_background;
			}

			if (strstr(mystring, "tukey.window.alpha") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atof(read_out.c_str()) != 0) {
						tukey_window_alpha = (float)atof(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ntukey.window.alpha = " << tukey_window_alpha << "\n";
			}

			if (strstr(mystring, "detector.width") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						detector_width = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ndetector.width = " << detector_width;
			}

			if (strstr(mystring, "detector.height") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						detector_height = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ndetector.height = " << detector_height << "\n";
			}

			if (strstr(mystring, "start.angle") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atof(read_out.c_str()) != 0) {
						start_angle = (float)atof(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nstart.angle = " << start_angle;
			}

			if (strstr(mystring, "end.angle") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atof(read_out.c_str()) != 0) {
						end_angle = (float)atof(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nend.angle = " << end_angle;
			}

			if (strstr(mystring, "projections.number") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						projections_number = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nprojections.number = " << projections_number << "\n";
			}

			if (strstr(mystring, "image.compression") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						image_compression = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nimage.compression = " << image_compression << "\n";
			}

			if (strstr(mystring, "voxels.per.cm") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						voxels_per_cm = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\nvoxels.per.cm = " << voxels_per_cm << "\n";
			}

			if (strstr(mystring, "gauss.blur.kernel") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						gauss_blur_kernel = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ngauss.blur.kernel = " << gauss_blur_kernel << "\n";
			}

			if (strstr(mystring, "a.multiplier") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atof(read_out.c_str()) != 0) {
						a_multiplier = (float)atof(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\na.multiplier = " << a_multiplier;
			}

			if (strstr(mystring, "gauss.noise.mean") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						gauss_noise_mean = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ngauss.noise.mean = " << gauss_noise_mean;
			}

			if (strstr(mystring, "gauss.noise.sigma") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						gauss_noise_sigma = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
				cout << "\ngauss.noise.sigma = " << gauss_noise_sigma << "\n";
			}

			if (strstr(mystring, "phantom.width") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						width = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "phantom.depth") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						depth = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "phantom.height") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						height = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "phantom.id") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						phantom_id = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "chamber.id") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						chamber_id = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}


			if (strstr(mystring, "save.phantom.yes") != NULL) {
				save_phantom_yes = 1;
			}

			if (strstr(mystring, "mammography.yes") != NULL) {
				mammography_yes = 1;
				cout << "\n*MAMMOGRAPHY MODE*\n";
			}

			if (strstr(mystring, "projections.yes") != NULL) {
				populate_projections = 1;
				cout << "\n*will generate the sequence of projections now*\n";
			}

			if (strstr(mystring, "defects.quantity") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						defects_quantity = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "defects.id") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						defects_id = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "voids.quantity") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						voids_quantity = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

			if (strstr(mystring, "voids.id") != NULL) {
				read_out = string(mystring);
				for (int i = 0; i < 200; ++i) {
					if (atoi(read_out.c_str()) != 0) {
						voids_id = atoi(read_out.c_str());
						break;
					}
					read_out.erase(0, 1);
				}
			}

		}
	}
	fclose(data_file);

	// отображение считанной информации для контроля правильности ввода в текстовом файле
	// displaying info on the phantom and random objects
	cout << "\n\nphantom.width=" << width << "\n";
	cout << "phantom.depth=" << depth << "\n";
	cout << "phantom.height=" << height << "\n";
	cout << "phantom.id=" << phantom_id << "\n";

	cout << "chamber.id=" << chamber_id << "\n";

	cout << "\ndefects.quantity=" << defects_quantity << "\n";
	cout << "defects.id=" << defects_id << "\n";

	cout << "\nvoids.quantity=" << voids_quantity << "\n";
	cout << "voids.id=" << voids_id << "\n";

	if (populate_projections == 1) {
		// Чтения файла спектра и нахождение количества ненулевых записей в нем
		// The following will create arrays to keep interpolated mu data and read the spectrum data
		GenerateMu muvalues(interpolation_points, spectrum_dir);

		// Создание фантома с заданными размерами
		// calling class Phantom constructor with phantom dimensions
		Phantom p_hantom(width, depth, height, chamber_id);
		// из выбранного материала
		// create an empty phantyom
		p_hantom.create_phantom(phantom_id);
		// выращивание случайных дефектов
		// growing the random defects
		p_hantom.grow_random_defects(defects_id, defects_quantity, voids_id, voids_quantity);

		// выращивание детерминированных дефектов
		// growing the custom made (deterministic) defects
		int object_presence = 0;
		int object_shape = 0, object_x = 0, object_y = 0, object_z = 0,
			object_width = 0, object_depth = 0, object_height = 0, element_id;
		float object_alpha = 0, object_beta = 0, object_gamma = 0;
		int object_id = 0;
		char shape_str[200], x_str[200], y_str[200], z_str[200],
			width_str[200], depth_str[200], height_str[200], alpha_str[200], beta_str[200], gamma_str[200],
			id_str[200], elementid[200];

		ofstream compare_mu_data("compare_mu_data.txt");
		cout << "\nParsed look-up table for element vs. id combination:\n\n";

		// продолжаем считывать файл input_file.txt
		// read data file again, now for the other parameters
		fopen_s(&data_file, "input_file.txt", "r");

		// подсчитаем количество элементов предпологаемых к использованию в фантоме
		// count the number of elements in the input_file.txt
		int elements_count = 0;

		// один достаточно большой, но ограниченный 256 элементами, массив для хранения данных по мю
		// keep the interpolated data for all elements (max=255) in a single variable
		int el_number = 255;
		float* interpolated_mudata = new float[el_number * (interpolation_points + 1)];
		float** interpolated_mu_data = new float*[el_number];
		// initialize it
		for (int i = 0; i < el_number; i++) {
			interpolated_mu_data[i] = interpolated_mudata + (interpolation_points + 1)*i;
		}

		if (data_file == NULL) perror("Error opening data file");
		else {
			while (!feof(data_file)) {
				fgets(mystring, 200, data_file);

				strcpy_s(shape_str, 200, mystring); strcpy_s(x_str, 200, mystring); strcpy_s(y_str, 200, mystring); strcpy_s(z_str, 200, mystring);
				strcpy_s(width_str, 200, mystring); strcpy_s(depth_str, 200, mystring); strcpy_s(height_str, 200, mystring);
				strcpy_s(alpha_str, 200, mystring); strcpy_s(beta_str, 200, mystring); strcpy_s(gamma_str, 200, mystring); strcpy_s(id_str, 200, mystring);

				strcpy_s(elementid, mystring);
				// считываем название элемента и номер поставленный ему в соответсвие
				// Parse look-up table for element vs. id combination
				if (strstr(mystring, "element.id") != NULL) {
					read_out = string(elementid);
					size_t found = read_out.find_first_of(":");
					int start_pos = int(found), end_pos = 0;
					string element_name;
					element_id = atoi(read_out.c_str());
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							element_id = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1); end_pos++;
					}
					if (element_id == 0) {
						element_name = "Vacuum";
					}
					else
					{
						read_out = string(mystring);
						read_out.erase(0, start_pos + 1);
						element_name = read_out.erase(end_pos - start_pos - 2, 200);

						muvalues.mu_by_name(element_name, d_ir);
						cout << element_name << "\t#" << element_id << "\n\n";

						// информация о энергии и соответствующем ему значении мю для данного элемента  содержится в 
						// muvalues.interp_mu[i][0] и muvalues.interp_mu[i][1] соответсвенно
						// 0-th row contains E values
						if (elements_count == 0) { // interpolation_points+1 to keep densities values later
							for (int i = 0; i < interpolation_points + 1; i++) {
								interpolated_mu_data[0][i] = muvalues.interp_mu[i][0];
							}
						}

						elements_count++; // interpolation_points+1 to keep densities values
						for (int i = 0; i < interpolation_points + 1; i++) {
							interpolated_mu_data[element_id][i] = muvalues.interp_mu[i][1];
						}

						// для просмотрщика скомпилированных спектров для используемых элементов в Матлаб
						// нужно удалить пробелы в названии элементов/веществ. Информация пишется в текстовой файл compare_mu_data.txt
						// need to remove empty spaces from the element's name. Needed for Matlab reader to build a legend..
						for (uchar i = 0; i < element_name.length(); i++) {
							if (element_name[i] == ' ') { element_name.erase(i, 1); }
						}

						compare_mu_data << element_name << "\t";
						for (int i = 0; i < interpolation_points; i++) {
							compare_mu_data << muvalues.interp_mu[i][0] << "\t";
						}
						compare_mu_data << "\n";
						compare_mu_data << element_name << "\t";
						for (int i = 0; i < interpolation_points; i++) {
							compare_mu_data << muvalues.interp_mu[i][1] << "\t";
						}
						compare_mu_data << "\n";
					}
				}
				// работаем с детерминированными объектами
				if (strstr(mystring, "object.shape") != NULL) {
					object_presence++;
					read_out = string(shape_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_shape = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.x") != NULL) {
					object_presence++;
					read_out = string(x_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_x = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.y") != NULL) {
					object_presence++;
					read_out = string(y_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_y = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.z") != NULL) {
					object_presence++;
					read_out = string(z_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_z = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.alpha") != NULL) {
					object_presence++;
					read_out = string(alpha_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_alpha = float(atoi(read_out.c_str())) / 180 * pi;
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.beta") != NULL) {
					object_presence++;
					read_out = string(beta_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_beta = float(atoi(read_out.c_str())) / 180 * pi;
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.gamma") != NULL) {
					object_presence++;
					read_out = string(gamma_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_gamma = float(atoi(read_out.c_str())) / 180 * pi;
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.width") != NULL) {
					object_presence++;
					read_out = string(width_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_width = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.depth") != NULL) {
					object_presence++;
					read_out = string(depth_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_depth = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.height") != NULL) {
					object_presence++;
					read_out = string(height_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_height = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}

				if (strstr(mystring, "object.id") != NULL) {
					object_presence++;
					read_out = string(id_str);
					for (int i = 0; i < 200; ++i) {
						if (atoi(read_out.c_str()) != 0) {
							object_id = atoi(read_out.c_str());
							break;
						}
						read_out.erase(0, 1);
					}
				}
				// описание КАЖДОГО детерминированного объекта должно иметь ТОЧНО 11 параметров
				if (object_presence == 11) { // object's description should have exactly 11 entries
					cout << "\nobject.shape = " << object_shape << "\n";
					cout << "object.x = " << object_x << "\n";
					cout << "object.y = " << object_y << "\n";
					cout << "object.z = " << object_z << "\n";
					cout << "object.width = " << object_width << "\n";
					cout << "object.depth = " << object_depth << "\n";
					cout << "object.height = " << object_height << "\n";
					cout << "object.alpha = " << object_alpha << "\n";
					cout << "object.beta = " << object_beta << "\n";
					cout << "object.gamma = " << object_gamma << "\n";
					cout << "object.id = " << object_id << "\n";

					p_hantom.o_bject(object_shape,
						object_x, object_y, object_z,
						object_width, object_depth, object_height,
						object_alpha, object_beta, object_gamma,
						object_id);

					object_presence = 0;
					object_shape = 0, object_x = 0, object_y = 0, object_z = 0,
						object_width = 0, object_depth = 0, object_height = 0;
					object_alpha = 0, object_beta = 0, object_gamma = 0;
					object_id = 0;
				}
			}
		}
		fclose(data_file);
		compare_mu_data.close();
		std::cout << "\n" << elements_count << " different elements/compounds are found in the table for processing\n";

		main_time = omp_get_wtime() - main_time;
		std::printf("\ntime to generate phantom in seconds = %f\n", main_time);

		if (save_phantom_yes == 1) {
			cout << "\nsaving the phantom's data ...";
			p_hantom.save_phantom();
		}
		else { cout << "\nno phantom saving" << "\n"; }

		// обрезаем часть пространство не содержащего фантом и не нуждающеося в пошаговой трассировке, см. ниже
		// attempt to optimize algorithm and cut the range of t values
		int max_dim = width;
		if (depth >= max_dim) { max_dim = depth; }
		if (height >= max_dim) { max_dim = height; }

		// гауссовское размытие и зашумление
		// set noise and blur values
		p_hantom.set_blur_and_noise(gauss_blur_kernel, a_multiplier, gauss_noise_mean, gauss_noise_sigma);

		/////////////////////////////////////////////////////
		cout << "\nproceeding to projections ...\n";

		// вычиляем размер шага который необходим для покрытия диапазона от начального угла поворота
		// до конечного, за столько то шагов равных числу проекций
		float degree_step;
		if (projections_number == 1) {
			degree_step = end_angle - start_angle;
		}
		else {
			degree_step = (end_angle - start_angle) / (float(projections_number) - 1.0f);
		}

		float t_before = 0, t_after = 1.0;
		int m_amo = 0;

		for (float d_egree = start_angle; d_egree <= end_angle; d_egree += degree_step) {
			if (mammography_yes == 0) {
				// обрезаем часть пространство не содержащего фантом и не нуждающеося в пошаговой трассировке, см. ниже
				int shift_degree = int(d_egree) % 90;
				t_before = (float(source_to_phantom) - (float)sin((45.0f + shift_degree) / 180.0f*pi)*float(max_dim) / 1.4142f) / float(source_to_phantom + phantom_to_detector);
				t_after = (float(source_to_phantom) + (float)sin((45.0f + shift_degree) / 180.0f*pi)*float(max_dim) / 1.4142f) / float(source_to_phantom + phantom_to_detector);

				p_hantom.create_projection(-source_to_phantom, source_elevation, phantom_to_detector, d_egree,
					detector_width, detector_height, interpolated_mu_data, muvalues.s_pectrum, interpolation_points,
					elements_count, muvalues.total_hw, image_compression, voxels_per_cm, t_before, t_after, m_amo, projections_dir);
				if ((end_angle == start_angle) | (projections_number == 1)) { break; }
			}
			else {
				m_amo = 1;
				// режим мамографии. Детектор вплотную придвинут к фантому и источник двигается по дуге окружности

				// phantom_to_detector = static_cast<int>(depth / 2.0f) + 1; // move one voxel further from the detector

				// в случае маммографии, старые переменные приобретают новый смысл и меняются с каждым шагом по углу
				// возвышения источника

				float source_to_detector_center = static_cast<float>(source_to_phantom + phantom_to_detector);
				int sourcetophantom = static_cast<int>(source_to_detector_center*cos(d_egree*pi / 180.0f)) - phantom_to_detector;
				int sourceelevation = static_cast<int>(source_to_detector_center*sin(d_egree*pi / 180.0f));

				t_before = (float(sourcetophantom) - float(depth) / 2.0f) / float(sourcetophantom + phantom_to_detector);
				t_after = (float(sourcetophantom) + float(depth) / 2.0f) / float(sourcetophantom + phantom_to_detector);

				p_hantom.create_projection(-sourcetophantom, sourceelevation, phantom_to_detector, d_egree,
					detector_width, detector_height, interpolated_mu_data, muvalues.s_pectrum, interpolation_points,
					elements_count, muvalues.total_hw, image_compression, voxels_per_cm, t_before, t_after, m_amo, projections_dir);

				if ((end_angle == start_angle) | (projections_number == 1)) { break; }
			}
		}


		delete[] interpolated_mu_data;
		delete[] interpolated_mudata;
	}

	//if(mammography_yes==1){
		// режим мамографии. Детектор вплотную придвинут к фантому и источник двигается по дуге окружности
		// с центром в центре детектора
	//	float source_to_detector_center = static_cast<float>(source_to_phantom + phantom_to_detector);
	//	int sourcetophantom = static_cast<int>(source_to_detector_center*cos(d_egree*pi / 180.0f) - phantom_to_detector);
	//	int sourceelevation = static_cast<int>(source_to_detector_center*sin(d_egree*pi / 180.0f));
	//}

		// процедура обратного проецирования проекций по заданым парметрам в input_file.txt
		cout << "\nproceeding to the supplied projections files ...\n";
		//for (int tcs = 190; tcs <= 200;tcs+=1) {
		//	target_cross_section=tcs;
			Loadprojection reconstruct_slice(-source_to_phantom, phantom_to_detector, detector_width, detector_height,
				start_angle, end_angle,
				F_OV, slice_plane, target_cross_section, slice_half_thickness, zero_padding_factor, cut_background, tukey_window_alpha,
				projections_number, projections_dir, mammography_yes);
			reconstruct_slice.load_images();
		//}

		std::cout << "done!\n";
	}
	// make some changes