//Do not edit the code below (unless you know what you are doing)
#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include"solution.h"
#include<fstream>
int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = matrix(L);
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(const matrix& A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double* A, int n)
{
	x = matrix(A, n);
	g = NAN;
	H = NAN;
	y = NAN;
}

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

ostream& operator<<(ostream& S, const solution& A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	S << "g_calls = " << solution::g_calls << endl;
	S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix O)
{
#if LAB_NO<2
	y = -1 * cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - 2 * 3.14, 2)) + 0.002 * (0.1 * x(0)) * (0.1 * x(0));
	++f_calls;
	matrix Y0(new double[3]{ 5,1,10 }, 3);
	matrix* Y = solve_ode(0, 1, 1000, Y0, x);



	ofstream file;
	file.open("FibY.txt");
	file << Y[1];
	file.close();

	double max = Y[1](0, 2);
	for (int i = 0; i < 1001; i++) {
		if (max < Y[1](i, 2)) {
			max = Y[1](i, 2);
		}

	}
#elif LAB_NO==3
	double a_ref = 3.14, o_ref = 0;
	matrix Y0(2, 1);
	matrix* Y = solve_ode(0, 0.1, 100, Y0, x); // x = [k1,k2]

	int* n = get_size(Y[1]);
	y(0) = 0;
	for (int i = 0; i < n[0]; ++i)
	{

		y(0) = y(0) + 10 * pow(a_ref - Y[1](i, 0), 2) + pow(o_ref - Y[1](i, 1), 2) + pow(x(0) * (a_ref - Y[1](i, 0)) + x(1) * (o_ref - Y[1](i, 1)), 2);
	}

	y(0) = y(0) * 0.1;
	/*y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;*/
	++f_calls;
#elif LAB_NO == 4
	matrix Y0(4, 1);
	Y0(0) = 0;
	Y0(1) = x(0);
	Y0(2) = 100;
	Y0(3) = 0;
	matrix* Y = solve_ode(0, 0.01, 7, Y0, x(1));
	double x0, x50;

	for (int i = 0; i < *get_size(Y[1]) - 1; i++) {
		if (fabs(Y[1](i, 2) - 50.0) < fabs(Y[1](i + 1, 2) - 50.0) && Y[1](i, 2) < 52 && Y[1](i, 2) > 48) {
			x50 = Y[1](i, 0);
		}
		if (Y[1](i, 2) > 0 && Y[1](i + 1, 2) < 0) {
			if (fabs(Y[1](i, 2)) < fabs(Y[1](i + 1, 2))) {
				x0 = Y[1](i, 0);
			}
			else {
				x0 = Y[1](i + 1, 0);
			}
		}
	}
	y = -x0;



	double g1, g2, g3, g4, g5, g6;
	g1 = -x(0) - 10;
	g2 = x(0) - 10;
	g3 = -x(1) - 20;
	g4 = x(1) - 20;
	g5 = -x50 + 4;
	g6 = x50 - 6;
	if (g1 > 0) {
		y = y + pow(g1, 2);
	}
	if (g2 > 0) {
		y = y + pow(g2, 2);
	}
	if (g3 > 0) {
		y = y + pow(g3, 2);
	}
	if (g4 > 0) {
		y = y + pow(g4, 2);
	}
	if (g5 > 0) {
		y = y + pow(g5, 2);
	}
	if (g6 > 0) {
		y = y + pow(g6, 2);
	}
	ofstream sym1("symulacja.txt");
	ofstream sym2("symulacja2.txt");
	for (int i = 0; i < *get_size(Y[1]); i++) {
		sym1 << Y[1](i, 0) << endl; //x
		sym2 << Y[1](i, 2) << endl; //y
	}
	++f_calls;
#elif LAB_NO == 5
	int* n = get_size(O);
	if (n[1] == 1) {
		y = pow((x(0) + 2 * x(1) - 7), 2) + pow(2 * x(0) + x(1) - 5, 2);
		++f_calls;
	}
	else {
		solution tmp;
		tmp.x = O[0] + x * O[1];
		tmp.fit_fun();
		y = tmp.y;
	}
#endif
}

#if LAB_NO==3
void solution::fit_fun_outside(matrix A)
{
	/*double arg = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
	y = sin(arg) / arg;

	if (-x(0) + 1 > 0)
		y = y + A(0) * pow(-x(0) + 1, 2);
	if (-x(1) + 1 > 0)
		y = y + A(0) * pow(-x(1) + 1, 2);
	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1) > 0)
		y = y + A(0) * pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1), 2);

	++f_calls;*/



	matrix Y0(4, 1);
	Y0(0) = 0;
	Y0(1) = x(0);
	Y0(2) = 100;
	Y0(3) = 0;
	matrix* Y = solve_ode(0, 0.01, 7, Y0, x(1));
	double x0, x50;

	for (int i = 0; i < *get_size(Y[1]) - 1; i++) {
		if (fabs(Y[1](i, 2) - 50.0) < fabs(Y[1](i + 1, 2) - 50.0) && Y[1](i, 2) < 52 && Y[1](i, 2) > 48) {
			x50 = Y[1](i, 0);
		}
		if (Y[1](i, 2) > 0 && Y[1](i + 1, 2) < 0) {
			if (fabs(Y[1](i, 2)) < fabs(Y[1](i + 1, 2))) {
				x0 = Y[1](i, 0);
			}
			else {
				x0 = Y[1](i + 1, 0);
			}
		}
	}
	y = -x0;



	double g1, g2, g3, g4, g5, g6;
	g1 = -x(0) - 10;
	g2 = x(0) - 10;
	g3 = -x(1) - 20;
	g4 = x(1) - 20;
	g5 = -x50 + 4;
	g6 = x50 - 6;
	if (g1 > 0) {
		y = y + pow(g1, 2);
	}
	if (g2 > 0) {
		y = y + pow(g2, 2);
	}
	if (g3 > 0) {
		y = y + pow(g3, 2);
	}
	if (g4 > 0) {
		y = y + pow(g4, 2);
	}
	if (g5 > 0) {
		y = y + pow(g5, 2);
	}
	if (g6 > 0) {
		y = y + pow(g6, 2);
	}
	ofstream sym("symulacja.txt");
	sym << "t x y" << endl;
	for (int i = 0; i < *get_size(Y[1]); i++) {
		sym << Y[1](i, 0) << " " << Y[1](i, 2) << endl;
	}

	++f_calls;
}

void solution::fit_fun_inside(matrix A)
{

	double arg = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
	y = sin(arg) / arg;

	if (-x(0) + 1 > 0)
		y = 1e10;
	else
		y = y - A(0) / (-x(0) + 1);
	if (-x(1) + 1 > 0)
		y = 1e10;
	else
		y = y - A(0) / (-x(1) + 1);
	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1) > 0)
		y = 1e10;
	else
		y = y - A(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1));

	++f_calls;
}
#endif

void solution::grad(matrix O)
{
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 10 * x(1) + 8 * x(0) - 38;
	++g_calls;
}

void solution::hess(matrix O)
{
	H = matrix(2, 2);
	H(0, 0) = 10;
	H(0, 1) = 8;
	H(1, 0) = 8;
	H(1, 1) = 10;
	++H_calls;
}
