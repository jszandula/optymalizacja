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
#endif
}

void solution::grad(matrix O)
{
	g = NAN;
	++g_calls;
}

void solution::hess(matrix O)
{
	H = NAN;
	++H_calls;
}
