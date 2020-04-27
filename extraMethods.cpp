#include "header.h"

vector<double> randomSearch(vector<double> x, int maxVal, int size, double delta0) {
	int count = 0;
	vector<vector<double>> steps(size, vector<double>(2));
	random_device random_device;
	mt19937 generator(random_device());

	do {
		double delta = delta0 - count / (maxVal / delta0);
		uniform_real_distribution<double> distribution1(x[0] - delta, x[0] + delta);
		uniform_real_distribution<double> distribution2(x[1] - delta, x[1] + delta);
		for (size_t i = 0; i < size; i++) {
			steps[i][0] = distribution1(generator);
			steps[i][1] = distribution2(generator);
		}
		x = argmin(argmin(steps, size), x);
		count++;
	} while (count < maxVal);

	return x;
}

vector<double> NelderMead(double eps) {
	double t = 5., s;
	vector<vector<double>> p(3, vector<double>(2));
	vector<double> r(2), r0(2);
	int xh, xg, xm;
	double d1 = t * (sqrt(3) + 1) / (2 * sqrt(2)), d2 = t * (sqrt(3) - 1) / (2 * sqrt(2));
	p[0][0] = 0;
	p[0][1] = 0;
	p[1][0] = d1;
	p[1][1] = d2;
	p[2][0] = d2;
	p[2][1] = d1;
	for (size_t i = 0; i < 3; i++) {
		r = sum(r, p[i]);
	}
	r = proizv(r, 1 / 3.);

	do {
		if (argmin(p[0], p[1], p[2]) == p[0]) {
			xm = 1;
			if (argmax(p[0], p[1], p[2]) == p[1]) {
				xh = 2;
				xg = 3;
			}
			else {
				xh = 3;
				xg = 2;
			}
		}
		else if (argmin(p[0], p[1], p[2]) == p[1]) {
			xm = 2;
			if (argmax(p[0], p[1], p[2]) == p[0]) {
				xh = 1;
				xg = 3;
			}
			else {
				xh = 3;
				xg = 1;
			}
		}
		else if (argmin(p[0], p[1], p[2]) == p[2]) {
			xm = 3;
			if (argmax(p[0], p[1], p[2]) == p[0]) {
				xh = 1;
				xg = 2;
			}
			else {
				xh = 2;
				xg = 1;
			}
		}

		vector<double> xc(2);
		for (size_t i = 0; i < 3; i++)
			if (i != (xh - 1))
				xc = sum(xc, p[i]);
		xc = proizv(xc, 0.5);

		vector<double> xr(2);
		xr = sum(proizv(xc, 2), proizv(p[xh - 1], -1));

		if (func(xr) < func(p[xm - 1])) {
			vector<double> xe(2);
			xe = sum(proizv(xc, -1), proizv(xr, 2));
			if (func(xe) < func(xr))
				p[xh - 1] = xe;
			else if (func(xe) > func(xr))
				p[xh - 1] = xr;
		}
		else if (func(xr) > func(p[xm - 1]) && func(xr) < func(p[xg - 1]))
			p[xh - 1] = xr;
		else if (func(xr) > func(p[xg - 1])) {
			if (func(xr) < func(p[xh - 1])) {
				vector<double> param(2);
				param = xr;
				xr = p[xh - 1];
				p[xh - 1] = param;
			}
			vector<double> xs(2);
			xs = sum(proizv(p[xh - 1], 0.5), proizv(xc, 0.5));
			if (func(xs) < func(p[xh - 1]))
				p[xh - 1] = xs;
			else if (func(xs) > func(p[xh - 1])) {
				for (size_t i = 0; i < 3; i++) {
					if (i != (xm - 1))
						p[i] = sum(p[xm - 1], proizv(sum(p[i], proizv(p[xm - 1], -1)), 0.5));
				}
			}
		}

		r0 = r;
		for (size_t i = 0; i < 3; i++) {
			r = sum(r, p[i]);
		}
		r = proizv(r, 1 / 3.);
		s = norm(sum(p[0], proizv(p[1], -1))) + norm(sum(p[0], proizv(p[2], -1))) + norm(sum(p[1], proizv(p[2], -1)));
	} while (norm(sum(r, proizv(r0, -1))) > eps && s > eps);

	return p[xm - 1];
}

vector<double> Powell(vector<double> x, double eps) {
	vector<double> x0(2);
	vector<vector<double>> p(2, vector<double>(2));

	p[0][0] = 1.;
	p[0][1] = 0.;
	p[1][0] = 0.;
	p[1][1] = 1.;

	do {
		x0 = x;
		for (size_t i = 0; i < 2; i++)
			x = squareInterp(x, p[i], 10, eps);
		vector<double> step(2);
		step = sum(proizv(x, 2), proizv(x0, -1));
		if (func(step) < func(x0)) {
			p[0] = p[1];
			p[1] = normalize(sum(x, proizv(x0, -1)));
			x = squareInterp(x, p[1], 10, eps);
		}
	} while (norm(sum(x, proizv(x0, -1))) > eps);

	return x;
}

vector<double> HookeJeeves(vector<double> x, double eps) {
	vector<double> h(2), x1(2), x2(2), x0(2);
	x1[0] = x[0];
	x1[1] = x[1];
	x2[0] = x[0];
	x2[1] = x[1];

	do {
		x0 = x;
		for (size_t i = 0; i < 2; i++) {
			h[i] = 4.;
			int count = 0;
			while (count == 0 && h[i] > eps) {
				x1[i] += h[i];
				x2[i] -= h[i];
				if (func(x1) < func(x) || func(x2) < func(x)) {
					x = argmin(x1, x2);
					count++;
				}
				else {
					x1[i] = x[i];
					x2[i] = x[i];
					h[i] /= 2.;
				}
			}
		}

		vector<double> p(2);
		for (size_t i = 0; i < 2; i++) {
			p[i] = x0[i] + 2 * (x[i] - x0[i]);
		}
		x = squareInterp(x, p, 10, eps);

	} while (h[0] > eps || h[1] > eps);

	return x;
}

vector<double> Rosenbrock(vector<double> x, double eps) {
	vector<double> x0, lambda(2);
	vector<vector<double>> p(2, vector<double>(2));

	p[0][0] = 1.;
	p[0][1] = 0.;
	p[1][0] = 0.;
	p[1][1] = 1.;

	do {
		x0 = x;

		for (size_t i = 0; i < 2; i++)
			lambda[i] = norm(sum(squareInterp(x, p[i], 10, eps), proizv(x, -1)));
		for (size_t i = 0; i < 2; i++)
			x = squareInterp(x, p[i], 10, eps);

		vector<vector<double>> a(2, vector<double>(2));
		a[1] = proizv(p[1], lambda[1]);
		a[0] = sum(a[1], proizv(p[0], lambda[0]));
		p[0] = normalize(a[0]);

		p[1] = normalize(sum(a[1], proizv(proizv(p[0], proizv(a[0], p[0])), -1)));
	} while (fabs(func(x0) - func(x)) > eps && fabs(norm(sum(x, proizv(x0, -1)))) > eps);

	return x;
}