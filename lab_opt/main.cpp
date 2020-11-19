#include<iostream>
#include<random>
#include<chrono>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"

using namespace std;

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==1
		double t0 = 0.0, dT = 0.1, theEnd = 50;
		matrix Y0 = matrix(new double[2]{ 1.0,2.0 }, 2);
		matrix* Y = solve_ode(t0, dT, theEnd, Y0);
		ofstream s("tout.csv");
		s << Y[0];
		s.close();
		s.open("yout.csv");
		s << Y[1];
		s.close();
#elif LAB_NO==2
		double x0, d, epsilion = 0.1, alfa;
		int Nmax = 60;
		double gamma = 0.01;
		d = 1;
		random_device R;
		ofstream exp_x0("exp_x0.txt");
		ofstream exp_a("exp_a.txt");
		ofstream exp_b("exp_b.txt");
		ofstream exp_calls("exp_calls.txt");

		ofstream fib_x("fib_x.txt");
		ofstream fib_y("fib_y.txt");
		ofstream fib_calls("fib_calls.txt");

		ofstream lag_x("lag_x.txt");
		ofstream lag_y("lag_y.txt");
		ofstream lag_calls("lag_calls.txt");

		//for (int i = 0; i < 100; i++) {
		x0 = 200.0 * R() / R.max() - 100;
		cout << "Pkt pocz.: " << x0 << endl;
		alfa = 2; // wsp.ekspacji - musi byæ wiêkszy od 1;
		exp_x0 << x0 << endl;
		cout << "Exp: ";
		double* p = expansion(x0, d, alfa, Nmax);

		cout << p[0] << " " << p[1];
		cout << endl;
		solution::clear_calls();
		exp_a << p[0] << endl;
		exp_b << p[1] << endl;


		// 0 - lokalne // 67 - globalne
		cout << "Fib: " << endl;
		solution x_fib = fib(-100, 100, epsilion);
		cout << x_fib << endl;
		fib_x << x_fib.x << endl;
		fib_y << x_fib.y << endl;
		fib_calls << x_fib.f_calls << endl;

		solution::clear_calls();

		cout << "Lag: " << endl;
		solution x_lag = lag(-100, 200, epsilion, gamma, Nmax);
		cout << x_lag;
		lag_x << x_lag.x << endl;
		lag_y << x_lag.y << endl;
		lag_calls << x_lag.f_calls << endl;

		solution::clear_calls;


		//}
		insert();
#elif LAB_NO==3
		/*solution test;
		test.x = matrix(new double[2]{ 1,1 }, 2);
		test.fit_fun();
		cout << test<<endl;*/
		ofstream x1_0("x1_0.txt");
		ofstream x2_0("x2_0.txt");

		ofstream HJ_x1("HJ_x1.txt");
		ofstream HJ_x2("HJ_x2.txt");
		ofstream HJ_y("HJ_y.txt");
		ofstream HJ_calls("HJ_calls.txt");

		ofstream Rosen_x1("Rosen_x1.txt");
		ofstream Rosen_x2("Rosen_x2.txt");
		ofstream Rosen_y("Rosen_y.txt");
		ofstream Rosen_calls("Rosen_calls.txt");

		matrix x0(2, 1);
		int Nmax = 1000;
		double s;			//poczatkowa dlugosc kroku
		double alfa;		//wpolczynnik skalowania dlugosci kroku
		double epsilon = 0.0001;	//dokladnosc
		double beta;
		matrix s0(2, 1);


		static default_random_engine generate(unsigned(time(nullptr)));
		uniform_real_distribution<double> distribution(-1, 1);

		//for (int i = 0; i < 100; i++) {
		x0(0) = 2.76993;//distribution(generate);
		x0(1) = 3.18721; //distribution(generate);
			std::cout << "Pkt pocz.: " << x0 << endl;
			x1_0 << x0(0) << endl;
			x2_0 << x0(1) << endl;

			s = 1.0;
			alfa = 0.5;
			std::cout << "HJ\n";
			solution x_HJ = HJ(x0, s, alfa, epsilon, Nmax);
			std::cout << x_HJ << endl;
			HJ_x1 << x_HJ.x(0) << endl;
			HJ_x2 << x_HJ.x(1) << endl;
			HJ_y << x_HJ.y << endl;
			HJ_calls << x_HJ.f_calls << endl;

			matrix Y0(2, 1);
			/*matrix* Y = solve_ode(0, 0.1, 100, Y0, x_HJ.x);
			ofstream S("Symulacja_HJ.txt");
			ofstream Sa("Symulacja_HJ_a.txt");
			for (int i = 0; i < 1000; i++) {
				Sa << Y[1][0](i) << endl;
				S << Y[1][1](i) << endl;

			}*/


			solution::clear_calls();

			s0(0) = s;
			s0(1) = s;
			alfa = 2.0;
			beta = 0.5;
			std::cout << "Rosenbrock\n";
			solution x_Rosen = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
			std::cout << x_Rosen;
			Rosen_x1 << x_Rosen.x(0) << endl;
			Rosen_x2 << x_Rosen.x(1) << endl;
			Rosen_y << x_Rosen.y << endl;
			Rosen_calls << x_Rosen.f_calls << endl;

			
			matrix* Y = solve_ode(0, 0.1, 100, Y0, x_Rosen.x);
			ofstream SR("Symulacja_R.txt");
			ofstream SRa("Symulacja_R_a.txt");
			for (int i = 0; i < 1000; i++) {
				SRa << Y[1][0](i) << endl;
				SR << Y[1][1](i) << endl;


			}

			solution::clear_calls;
		//}

	
#endif
	}
	catch (char* EX_INFO)
	{
		std::cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
