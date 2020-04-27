#include <iostream>
#include <vector>
#include <random>
#include <fstream>

using namespace std;

double func(vector<double> x);
double norm(vector<double> vec);
vector<double> normalize(vector<double> vec);
vector<double> sum(vector<double> vec1, vector<double> vec2);
vector<double> proizv(vector<double> vec, double num);
double proizv(vector<double> vec1, vector<double> vec2);
vector<double> argmin(vector<double> x1, vector<double> x2, vector<double> x3);
vector<double> argmin(vector<double> x1, vector<double> x2);
vector<double> argmax(vector<double> x1, vector<double> x2, vector<double> x3);
vector<double> argmin(vector<vector<double>> p, int size);

vector<double> randomSearch(vector<double> x, int maxVal = 1000, int size = 6, double delta0 = 10);
vector<double> NelderMead(double eps);
vector<double> Powell(vector<double> x, double eps);
vector<double> HookeJeeves(vector<double> x, double eps);
vector<double> Rosenbrock(vector<double> x, double eps);

vector<double> slidingWindow(vector<double> x, vector<double> p, double h);
vector<double> squareInterp(vector<double> x, vector<double> p, double h0, double eps);
vector<double> coordinateDescent(vector<double> x, double eps);