//Do not edit the code below (unless you know what you are doing)

#define _USE_MATH_DEFINES
#include<cmath>
#include"ode_solver.h"


matrix* solve_ode(double t0, double dt, double tend, const matrix& Y0, matrix P)
{
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	if (N < 2)
		throw "The time interval is defined incorrectly";
	int* s = get_size(Y0);
	if (s[1] != 1)
		throw "Initial condition must be a vector";
	int n = s[0];
	delete[]s;
	matrix* S = new matrix[2]{ matrix(N), matrix(n,N) };
	S[0](0) = t0;
	for (int i = 0; i < n; ++i)
		S[1](i, 0) = Y0(i);
	matrix k1(n), k2(n), k3(n), k4(n);
	for (int i = 1; i < N; ++i)
	{
		S[0](i) = S[0](i - 1) + dt;
		k1 = dt * diff(S[0](i - 1), S[1][i - 1], P);
		k2 = dt * diff(S[0](i - 1) + 0.5 * dt, S[1][i - 1] + 0.5 * k1, P);
		k3 = dt * diff(S[0](i - 1) + 0.5 * dt, S[1][i - 1] + 0.5 * k2, P);
		k4 = dt * diff(S[0](i - 1) + dt, S[1][i - 1] + k3, P);
		for (int j = 0; j < n; ++j)
			S[1](j, i) = S[1](j, i - 1) + (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j)) / 6;
	}
	S[1] = trans(S[1]);
	return S;
}

//You can edit the following code

matrix diff(double t, const matrix& Y, matrix P)
{
#if LAB_NO==1
	double m = 5, b = 1.5, k = 1, f = 0;
	matrix dY(Y);
	dY(0) = Y(1);
	dY(1) = (f - b * Y(1)) - k * Y(0)) / m;
	return dY;
#elif LAB_NO == 2
	//matrix dY;
	////JEDNOSTKI np Fin i DB 
	//double a = 0.98, b = 0.63, g = 9.81;
	//double PA = 1, PB = 1, DB = 36.5665;
	//double Fin = 10, Tin = 10, TA = 90;
	//double DA = P(0);	//x(0)
	//double FAout = Y(0) > 0 ? a * b * DA * sqrt(2 * g * Y(0) / PA) : 0;
	//double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	//dY(0) = -FAout; //woda wylewajaca sie ze zbiornika A
	//dY(1) = FAout + Fin - FBout; //woda wlewajaca sie do B
	////dY(2) = Vin/V (Tin - T) jak zmienia sie temperatura
	//return dY;

	matrix dY(Y);

	double a, b, g, PA, PB, DB, Fin, Tin, TA, TB, DA, FAout, FBout, TBout, TBin;
	a = 0.98;
	b = 0.63;
	g = 9.81;
	PA = 1;
	TA = 90;
	TB = 10;
	PB = 1;
	Tin = 10;
	Fin = 0.001; //10 l/s
	DB = 0.0003665;//36.5665 cm2;


	DA = P(0); // P to bedzie wielkosci otworu p to maciesz 1x1 zawierajaca wartosc zmiennej decyzyjnej DA

	double VAout = Y(0) <= 0 ? 0.0 : -a * b * DA * sqrt(2 * g * Y(0) / PA);
	//ilo�� wody wyp�ywaj�ca ze zbiornika A
	double VBout = Y(1) <= 0 ? 0.0 : -a * b * DB * sqrt(2 * g * Y(1) / PB);

	// TA 90
	TBin = (((-VAout * 1000) * TA) + (10 * 10)) / ((-VAout * 1000) + 10);
	TBout = ((Fin - VAout) / Y(1)) * (TBin - Y(2));

	dY(0) = VAout;		//mozliwe ze bez minusa 
	dY(1) = -VAout + Fin + VBout;

	dY(2) = TBout;

	return dY;


	// to nie dziala
	/*
	dY(0) = -(Y(0) > 0 ? a * b * DA * sqrt(2 * g * Y(0) / PA) : 0);
	dY(1) = -dY(0) + Fin - (Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0);
	dY(2) = Fin / Y(1) * (Tin - Y(2)) - dY(0) / Y(1) * (TA - Y(2));*/
#elif LAB_NO==3
	double mr = 1, mc = 10, L = 0.5, b = 0.5, a_ref, o_ref;
	a_ref = M_PI;
	o_ref = 0;

	double I = mc * L * L + mr * L * L / 3;
	double k1 = P(0), k2 = P(1);
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (k1 * (a_ref - Y(0)) + k2 * (o_ref - Y(1)) - b * Y(1)) / I;

	return dY;
#elif LAB_NO == 4
	double C = 0.47, r = 0.12, m = 0.6, ro = 1.2, g = 9.81;
	double S = 3.1415 * r * r;
	double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));
	//p(0) - omega
	double Fmx = 3.1415 * ro * Y(3) * P(0) * pow(r, 3);
	double Fmy = 3.1415 * ro * Y(1) * P(0) * pow(r, 3);
	matrix dY(Y);
	dY(0) = Y(1);
	dY(1) = (-Dx - Fmx) / m;
	dY(2) = Y(3);
	dY(3) = (-(m * g) - Dy - Fmy) / m;

	return dY;
#elif LAB_NO ==7
	double m1 = 5.0, m2 = 5.0, k1 = 1.0, k2 = 1.0, F = 1.0;
	double b1 = P(0, 0), b2 = P(1, 0);
	matrix dY(4, 1);
	dY(0) = Y(1);
	dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
	dY(2) = Y(3);
	dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;
	return dY;


#else
	matrix dY;
	return dY;
#endif
}