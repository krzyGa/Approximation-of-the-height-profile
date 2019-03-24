#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <time.h>

#define SIZE_ORG 500


using namespace std;

struct value {
	double x;
	double f;
};

value* load_data(string path) {
	string line;
	ifstream file;
	file.open(path);

	if (!file.is_open()) {
		cout << "Plik niedostepny";
		return nullptr;
	}
	value *xf = new value[SIZE_ORG];

	double d, e;
	char c;
	for (int i = 0; i < SIZE_ORG; i++) {
		if (!getline(file, line)) break;

		file >> d;
		file >> c;
		file >> e;

		xf[i].x = d;
		xf[i].f = e;
	}

	file.close();
	return xf;
}

value* mofidy_loaded_data(value *v, int number_of_points, bool evenly)
{
	int space = SIZE_ORG / number_of_points;
	int size = number_of_points;
	value *v_return = new value[size];
	int j = 0;
	if (evenly) {
		for (int i = 0; i < SIZE_ORG; i++) {
			if (i % space == 0) {
				v_return[j].x = v[i].x;
				v_return[j].f = v[i].f;
				j++;
			}
		}
		--j;
		v_return[j].x = v[SIZE_ORG - 1].x;
		v_return[j].f = v[SIZE_ORG - 1].f;
	}
	else {
		int space_ran_max = space;
		v_return[0].x = v[0].x;
		v_return[0].f = v[0].f;
		for (int i = 1; i < number_of_points-1; i++) {
			int temp = rand() % (space_ran_max + 1 - space_ran_max - space) + (space_ran_max - space);
			cout << space_ran_max << " " << (space_ran_max - space) << endl;
			v_return[i].x = v[temp].x;
			v_return[i].f = v[temp].f;
			space_ran_max += space;
		}
		v_return[number_of_points - 1].x = v[SIZE_ORG - 1].x;
		v_return[number_of_points - 1].f = v[SIZE_ORG - 1].f;
	}


	return v_return;
}

void save_data(string fileName, string colName, value *v, int size) {
	ofstream toFile;
	toFile.open(fileName);

	toFile << "X " << colName << std::endl;
	
	for (int i = 0; i < size; i++) {
		toFile << v[i].x << " " << v[i].f << endl;
	}

	toFile.close();
}

void thomas_method(double *a, double *b, double *c, double *x, double *d, unsigned int n) {
	for (int i = 1; i < n; i++) {
		double l = a[i] / b[i - 1];
		a[i] = 0.0;
		b[i] -= l * c[i - 1];
		d[i] -= l * d[i - 1];
	}

	x[n - 1] = d[n - 1] / b[n - 1];
	for (int i = n - 2; i >= 0; i--) x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
}

double spline_interpolation(const value *xf, unsigned int xf_size, double x) {
	double *_a = new double[xf_size];
	double *_b = new double[xf_size];
	double *_c = new double[xf_size];
	double *_x = new double[xf_size];
	double *_d = new double[xf_size];
	int i0 = 0;
	bool f = false;

	for (int i = 1; i < xf_size - 1; i++) {
		double dxi = xf[i + 1].x - xf[i].x;
		double dxi1 = xf[i].x - xf[i - 1].x;
		_a[i] = dxi1 / dxi;
		_b[i] = 2 * (xf[i + 1].x - xf[i - 1].x) / (dxi);
		_c[i] = 1;
		_d[i] = 6 * (((xf[i + 1].f - xf[i].f) / (dxi * dxi)) - ((xf[i].f
			- xf[i - 1].f) / (dxi * dxi1)));

		if (!f && xf[i].x <= x && x <= xf[i + 1].x) {
			i0 = i;
			f = true;
		}
	}

	_x[0] = _x[xf_size - 1] = 0;
	thomas_method(_a + 1, _b + 1, _c + 1, _x + 1, _d + 1, xf_size - 2);

	double dxi = xf[i0 + 1].x - xf[i0].x;
	double dxfi0 = x - xf[i0].x;
	double dxfi1 = xf[i0 + 1].x - x;
	double fx = (_x[i0] / 6.0) * ((dxfi1 * dxfi1 * dxfi1) / dxi - dxi * dxfi1)
		+ (_x[i0 + 1] / 6.0) * ((dxfi0 * dxfi0 * dxfi0) / dxi - dxi * dxfi0)
		+ xf[i0].f * dxfi1 / dxi + xf[i0 + 1].f * dxfi0 / dxi;

	delete[] _d;
	delete[] _x;
	delete[] _c;
	delete[] _b;
	delete[] _a;
	return fx;
}

void save_spline_interpolation_data(string fileName, string colName, value *v, int size) {
	ofstream toFile;
	toFile.open(fileName);

	toFile << "X " << colName << std::endl;

	for (double x = v[0].x; x < v[size-1].x; x += 1.0) {
		toFile << x << " " << spline_interpolation(v, size, x) << std::endl;
	}

	toFile.close();
}

double lagrange_interpolation(const value *xf, unsigned int xf_size, double x) {
	double result = 0.0;
	for (int i = 0; i < xf_size; i++) {
		double f1 = 1.0;
		double f2 = 1.0;
		for (int j = 0; j < xf_size; j++){
			if (i != j) {
				f1 *= x - xf[j].x;
				f2 *= xf[i].x - xf[j].x;
				}
		}
		result += xf[i].f * f1 / f2;
	}
	return result;
}

void save_lagrange_interpolation_data(string fileName, string colName, value *v, int size) {
	ofstream toFile;
	toFile.open(fileName);

	toFile << "X " << colName << std::endl;

	for (double x = v[0].x; x < v[size - 1].x; x += 1.0) {
		toFile << x << " " << lagrange_interpolation(v, size, x) << std::endl;
	}

	toFile.close();
}


int main()
{
	srand(time(NULL));

	// in order to generete out for another file just change its name in load_data
	value *xf_org = load_data("no_hill.csv");
	save_data("data_org.dat", "$ORG$", xf_org, SIZE_ORG);

	int n = 10;
	value *xf_mod_evenly = mofidy_loaded_data(xf_org, n, true);
	save_data("data_mod_n_10.dat", "$MOD$", xf_mod_evenly, n);
	save_spline_interpolation_data("data_n_10.dat", "$N_10$", xf_mod_evenly, n);
	save_lagrange_interpolation_data("data_nL_10.dat", "$N_10$", xf_mod_evenly, n);


	n = 5;
	xf_mod_evenly = mofidy_loaded_data(xf_org, n, true);
	//save_spline_interpolation_data("data_n_5.dat", "$N_5$", xf_mod_evenly, n);
	save_data("data_mod_n_5.dat", "$MOD$", xf_mod_evenly, n);
	save_lagrange_interpolation_data("data_nL_5.dat", "$N_5$", xf_mod_evenly, n);

	n = 20;
	xf_mod_evenly = mofidy_loaded_data(xf_org, n, true);
	//save_spline_interpolation_data("data_n_20.dat", "$N_20$", xf_mod_evenly, n);
	save_lagrange_interpolation_data("data_nL_20.dat", "$N_20$", xf_mod_evenly, n);


	n = 25;
	xf_mod_evenly = mofidy_loaded_data(xf_org, n, true);
	//save_data("data_mod_n_25.dat", "$MOD$", xf_mod_evenly, n);
	save_spline_interpolation_data("data_n_25.dat", "$N_25$", xf_mod_evenly, n);
	//save_lagrange_interpolation_data("data_nL_25.dat", "$N_25$", xf_mod_evenly, n);


	n = 50;
	xf_mod_evenly = mofidy_loaded_data(xf_org, n, true);
	//save_data("data_mod_n_50.dat", "$MOD$", xf_mod_evenly, n);
	save_spline_interpolation_data("data_n_50.dat", "$N_50$", xf_mod_evenly, n);
	//save_lagrange_interpolation_data("data_nL_50.dat", "$N_50$", xf_mod_evenly, n);


	n = 10;
	value *xf_mod_random = mofidy_loaded_data(xf_org, n, false);
	save_data("data_mod_ran_n_10.dat", "$RAN$", xf_mod_random, n);
	//save_spline_interpolation_data("data_ran_n_10.dat", "$N_10$", xf_mod_random, n);
	save_lagrange_interpolation_data("data_ran_nL_10.dat", "$N_10$", xf_mod_random, n);

	n = 5;
	xf_mod_random = mofidy_loaded_data(xf_org, n, false);
	save_data("data_mod_ran_n_5.dat", "$RAN$", xf_mod_random, n);
	//save_spline_interpolation_data("data_ran_n_5.dat", "$N_5$", xf_mod_random, n);
	save_lagrange_interpolation_data("data_ran_nL_5.dat", "$N_5$", xf_mod_random, n);


	n = 25;
	xf_mod_random = mofidy_loaded_data(xf_org, n, false);
	//save_data("data_mod_ran_n_25.dat", "$RAN$", xf_mod_random, n);
	save_spline_interpolation_data("data_ran_n_25.dat", "$N_25$", xf_mod_random, n);
	//save_lagrange_interpolation_data("data_ran_nL_25.dat", "$N_25", xf_mod_random, n);


	n = 50;
	xf_mod_random = mofidy_loaded_data(xf_org, n, false);
	//save_data("data_mod_ran_n_50.dat", "$RAN$", xf_mod_random, n);
	save_spline_interpolation_data("data_ran_n_50.dat", "$N_50$", xf_mod_random, n);
	//save_lagrange_interpolation_data("data_ran_nL_50.dat", "$N_50", xf_mod_random, n);

	return 0;
}