#include "header.h"

vector<double> slidingWindow(vector<double> x, vector<double> p, double h) {
	vector<double> x1(2), x2(2), interval(2);

	if (p[0] != 0) {
		for (;;) {
			x1[0] = x[0];
			x1[1] = x[1];
			x2[0] = x[0];
			x2[1] = x[1];
			x1[0] -= h;
			x2[0] += h;
			x1[1] -= p[1] / p[0] * h;
			x2[1] += p[1] / p[0] * h;
			interval[0] = x1[0];
			interval[1] = x2[0];

			if (func(x1) > func(x) && func(x) < func(x2))
				return interval;
			else if (func(x1) > func(x2)) {
				x[0] += h / 2.;
				x[1] += p[1] / p[0] * h / 2.;
			}
			else if (func(x1) < func(x2)) {
				x[0] -= h / 2.;
				x[1] -= p[1] / p[0] * h / 2.;
			}
		}
	}
	else {
		for (;;) {
			x1[0] = x[0];
			x1[1] = x[1];
			x2[0] = x[0];
			x2[1] = x[1];
			x1[1] -= h;
			x2[1] += h;
			interval[0] = x1[1];
			interval[1] = x2[1];

			if (func(x1) > func(x) && func(x) < func(x2))
				return interval;
			else if (func(x1) > func(x2))
				x[1] += h / 2.;
			else if (func(x1) < func(x2))
				x[1] -= h / 2.;
		}
	}
}

vector<double> squareInterp(vector<double> x, vector<double> p, double h0, double eps) {
	double h = 2.;
	vector<double> x1(2), x2(2), x3(2), X(2);
	x1[0] = x[0];
	x1[1] = x[1];
	x2[0] = x[0];
	x2[1] = x[1];
	x3[0] = x[0];
	x3[1] = x[1];
	X[0] = x[0];
	X[1] = x[1];
	vector<double> interval(2);
	interval = slidingWindow(x, p, h0);
	vector<double> xmin(2);


	if (p[0] != 0) {
		x1[0] = (interval[0] + interval[1]) / 2;
		x1[1] += p[1] / p[0] * (x1[0] - x[0]);
		do {
			x2[0] = x1[0] + h;
			x2[1] += p[1] / p[0] * (x2[0] - x[0]);
			if (func(x1) > func(x2))
				x3[0] = x1[0] + 2 * h;
			else
				x3[0] = x1[0] - h;
			x3[1] += p[1] / p[0] * (x3[0] - x[0]);
			X[0] = (x2[0] * x2[0] - x3[0] * x3[0]) * func(x1) + (x3[0] * x3[0] - x1[0] * x1[0]) * func(x2) + (x1[0] * x1[0] - x2[0] * x2[0]) * func(x3);
			X[0] /= 2 * ((x2[0] - x3[0]) * func(x1) + (x3[0] - x1[0]) * func(x2) + (x1[0] - x2[0]) * func(x3));
			X[1] += p[1] / p[0] * (X[0] - x[0]);
			xmin = argmin(x1, x2, x3);
			if (X[0] >= x1[0] && X[0] <= x3[0])
				x1 = argmin(xmin, X);
			else
				x1[0] = X[0];
			x1[1] += p[1] / p[0] * (x1[0] - x[0]);
		} while (fabs(func(xmin) - func(X)) > eps && fabs(xmin[0] - X[0]) > eps);
	}
	else {
		x1[1] = (interval[0] + interval[1]) / 2;
		do {
			x2[1] = x1[1] + h;
			if (func(x1) > func(x2))
				x3[1] = x1[1] + 2 * h;
			else
				x3[1] = x1[1] - h;
			X[1] = (x2[1] * x2[1] - x3[1] * x3[1]) * func(x1) + (x3[1] * x3[1] - x1[1] * x1[1]) * func(x2) + (x1[1] * x1[1] - x2[1] * x2[1]) * func(x3);
			X[1] /= 2 * ((x2[1] - x3[1]) * func(x1) + (x3[1] - x1[1]) * func(x2) + (x1[1] - x2[1]) * func(x3));
			xmin = argmin(x1, x2, x3);
			if (X[1] >= x1[1] && X[1] <= x3[1])
				x1 = argmin(xmin, X);
			else
				x1[1] = X[1];
		} while (fabs(func(xmin) - func(X)) > eps && fabs(xmin[1] - X[1]) > eps);
	}

	return X;
}

vector<double> coordinateDescent(vector<double> x, double eps) {
	int iter = 0, coef;
	vector<double> x0, p(2);
	ofstream ofs;
	ofs.open("coordinateDescent.txt");
	ofs << x[0] << " " << x[1] << endl;
	do {
		x0 = x;
		coef = (iter % 2 == 0) ? 0 : 1;
		p[coef] = 1;
		x = squareInterp(x, p, 10, eps);
		iter++;
		p[coef] = 0;
		ofs << x[0] << " " << x[1] << endl;
	} while (fabs(func(x0) - func(x)) > eps && fabs(x0[coef] - x[coef]) > eps);
	ofs.close();
	cout << iter << " iterations\n";

	return x;
}

int main(void) {
	const double eps = 1e-7;
	vector<double> x(2);
	x[0] = 9.;
	x[1] = 25.;
	vector<double> extr(2);

	cout << "Coordinate descent method\n";
	extr = coordinateDescent(x, eps);
	cout << "Maximum is on point (" << extr[0] << ", " << extr[1] << ")\nfmax = " << -func(extr) << endl << endl;

	extr = randomSearch(x);
	cout << "Random search method (" << extr[0] << ", " << extr[1] << ")\n"; 
	
	extr = NelderMead(eps);
	cout << "Nelder - Mead method (" << extr[0] << ", " << extr[1] << ")\n";

	extr = Powell(x, eps);
	cout << "Powell method (" << extr[0] << ", " << extr[1] << ")\n";

	extr = HookeJeeves(x, eps);
	cout << "Hooke - Jeeves method (" << extr[0] << ", " << extr[1] << ")\n";

	extr = Rosenbrock(x, eps);
	cout << "Rosenbrock method (" << extr[0] << ", " << extr[1] << ")\n";

	return 0;
}
