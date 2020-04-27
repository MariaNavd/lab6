#include "header.h"

double func(vector<double> x) {
	return x[0] * x[0] / 3. + x[0] * x[1] / 6. - 127 * x[0] / 6. + x[1] * x[1] / 3. - 91 * x[1] / 6. + 239 / 2.;
}

double norm(vector<double> vec) {
	double sum = 0;
	for (size_t i = 0; i < vec.size(); i++)
		sum += vec[i] * vec[i];
	return sqrt(sum);
}

vector<double> normalize(vector<double> vec) {
	for (size_t i = 0; i < vec.size(); i++)
		vec[i] /= norm(vec);
	return vec;
}

vector<double> sum(vector<double> vec1, vector<double> vec2) {
	vector<double> sum;
	for (size_t i = 0; i < vec1.size(); i++) {
		sum.push_back(vec1[i] + vec2[i]);
	}
	return sum;
}

vector<double> proizv(vector<double> vec, double num) {
	vector<double> comp;
	for (size_t i = 0; i < vec.size(); i++) {
		comp.push_back(num * vec[i]);
	}
	return comp;
}

double proizv(vector<double> vec1, vector<double> vec2) {
	double sum = 0;
	for (size_t i = 0; i < vec1.size(); i++) {
		sum += vec1[i] * vec2[i];
	}
	return sum;
}

vector<double> argmin(vector<double> x1, vector<double> x2, vector<double> x3) {
	if (func(x1) < func(x2)) {
		if (func(x1) < func(x3))
			return x1;
		else
			return x3;
	}
	else {
		if (func(x2) < func(x3))
			return x2;
		else
			return x3;
	}
}

vector<double> argmin(vector<double> x1, vector<double> x2) {
	if (func(x1) < func(x2))
		return x1;
	else
		return x2;
}

vector<double> argmax(vector<double> x1, vector<double> x2, vector<double> x3) {
	if (func(x1) > func(x2)) {
		if (func(x1) > func(x3))
			return x1;
		else
			return x3;
	}
	else {
		if (func(x2) > func(x3))
			return x2;
		else
			return x3;
	}
}

vector<double> argmin(vector<vector<double>> p, int size) {
	vector<double> mVal;
	mVal = p[0];
	for (size_t i = 0; i < size; i++) {
		if (func(p[i]) < func(mVal))
			mVal = p[i];
	}
	return mVal;
}