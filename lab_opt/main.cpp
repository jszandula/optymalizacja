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

#elif LAB_NO==4
		ofstream f("lab4_results.csv");
		
		ofstream r_out_file("niewiem1.txt");
		ofstream r_in_file("niewiem2.txt");
		
		/*
		double epsilon = 1e-5, c = 1000, a;
		int Nmax = 10000;
		matrix x0(2, 1);
		random_device rd;
		a = 4;*/
		/*
		for (size_t i = 0; i < 100; i++)
		{
			cout << "a = " << a << " ; iteration: " << i << endl;
			do
			{
				x0(0) = 4.0*rd() / rd.max() + 1;
				x0(1) = 4.0*rd() / rd.max() + 1;
			} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);
			solution opt_out = pen_outside(x0, c, a, epsilon, Nmax);
			double r_out = sqrt(pow(opt_out.x(0), 2) + pow(opt_out.x(1), 2));
			double f_calls_out = opt_out.f_calls;
			solution::clear_calls();
			solution opt_in = pen_inside(x0, c, a, epsilon, Nmax);
			double r_in = sqrt(pow(opt_in.x(0), 2) + pow(opt_in.x(1), 2));
			double f_calls_in = opt_in.f_calls;
			solution::clear_calls();
			f << x0(0) << ";" << x0(1) << ";" << opt_out.x(0) << ";" << opt_out.x(1) << ";" << r_out << ";" << opt_out.y(0) << ";" << f_calls_out
			 << ";" << opt_in.x(0) << ";" << opt_in.x(1) << ";" << r_in << ";" << opt_in.y(0) << ";" << f_calls_in << ";" << endl;
		}
		a = 4.4934;
		for (size_t i = 0; i < 100; i++)
		{
			cout << "a = " << a << " ; iteration: " << i << endl;
			do
			{
				x0(0) = 4.0*rd() / rd.max() + 1;
				x0(1) = 4.0*rd() / rd.max() + 1;
			} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);
			solution opt_out = pen_outside(x0, c, a, epsilon, Nmax);
			double r_out = sqrt(pow(opt_out.x(0), 2) + pow(opt_out.x(1), 2));
			r_out_file << r_out << endl;
			double f_calls_out = opt_out.f_calls;
			solution::clear_calls();
			solution opt_in = pen_inside(x0, c, a, epsilon, Nmax);
			double r_in = sqrt(pow(opt_in.x(0), 2) + pow(opt_in.x(1), 2));
			r_in_file << r_in << endl;
			double f_calls_in =opt_in.f_calls;
			solution::clear_calls();
			f << x0(0) << ";" << x0(1) << ";" << opt_out.x(0) << ";" << opt_out.x(1) << ";" << r_out << ";" << opt_out.y(0) << ";" << f_calls_out
			  << ";" << opt_in.x(0) << ";" << opt_in.x(1) << ";" << r_in << ";" << opt_in.y(0) << ";" << f_calls_in << ";" << endl;
		}
		a = 5;
		for (size_t i = 0; i < 100; i++)
		{
			cout << "a = " << a << " ; iteration: " << i << endl;
			do
			{
				x0(0) = 4.0*rd() / rd.max() + 1;
				x0(1) = 4.0*rd() / rd.max() + 1;
			} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);
			solution opt_out = pen_outside(x0, c, a, epsilon, Nmax);
			double r_out = sqrt(pow(opt_out.x(0), 2) + pow(opt_out.x(1), 2));
			r_out_file << r_out << endl;
			double f_calls_out = opt_out.f_calls;
			solution::clear_calls();
			solution opt_in = pen_inside(x0, c, a, epsilon, Nmax);
			double r_in = sqrt(pow(opt_in.x(0), 2) + pow(opt_in.x(1), 2));
			r_in_file << r_in << endl;
			double f_calls_in = opt_in.f_calls;
			solution::clear_calls();
			f << x0(0) << ";" << x0(1) << ";" << opt_out.x(0) << ";" << opt_out.x(1) << ";" << r_out << ";" << opt_out.y(0) << ";" << f_calls_out
			  << ";" << opt_in.x(0) << ";" << opt_in.x(1) << ";" << r_in << ";" << opt_in.y(0) << ";" << f_calls_in << ";" << endl;
		}*/

		matrix x0(2, 1);
		x0(0) = 2.1;
		x0(1) = 6;
		double c0 = 1.0, dc = 2.0, epsilon = 0.0001;
		int Nmax = 5000;
		matrix x1(2, 1);
		solution sympleks = pen(x0, c0, dc, epsilon, Nmax, x1);
		
		std::cout << "Optymalne V0x= " << sympleks.x(0) << "\nOptymalna omega: " << sympleks.x(1);
		std::cout << endl << "x koncowe: " << (-1) * sympleks.y(0) << endl;
		std::cout << "f_calls: " << sympleks.f_calls << endl;

		solution::clear_calls();
#elif LAB_NO == 5
/*
		ofstream f("lab5_res.csv");
		double epsilon = 0.0001;
		int Nmax = 1000, SDf_calls, SDg_calls, CGf_calls, CGg_calls, Newf_calls, Newg_calls, Newh_calls,
					SDf_calls1, SDg_calls1, CGf_calls1, CGg_calls1, Newf_calls1, Newg_calls1, Newh_calls1, 
					SDf_calls2, SDg_calls2, CGf_calls2, CGg_calls2, Newf_calls2, Newg_calls2, Newh_calls2;
		matrix limits(2, 2);
		limits(0,0) = -10;
		limits(0,1) = 10;
		limits(1,0) = -10;
		limits(1,1) = 10;
		matrix x0(2, 1);

		
		x0(0) = 7.24004;
		x0(1) = -1.72375;


		cout << "\nSD 0,05" << endl;
		//h stalokrokowe 0,05
		solution SDres = SD(x0, 0.05, epsilon, Nmax, limits);
		SDf_calls = solution::f_calls;
		SDg_calls = solution::g_calls;
		solution::clear_calls();

		cout << "CG 0,05" << endl;
		//h stalokrokowe 0,05
		solution CGres = CG(x0, 0.05, epsilon, Nmax, limits);
		CGf_calls = solution::f_calls;
		CGg_calls = solution::g_calls;
		solution::clear_calls();

		cout << "Newton 0,05" << endl;
		//h stalokrokowe 0,05
		solution Newtonres = Newton(x0, 0.05, epsilon, Nmax, limits);
		Newf_calls = solution::f_calls;
		Newg_calls = solution::g_calls;
		Newh_calls = solution::H_calls;
		solution::clear_calls();

		cout << "SD 0,12" << endl;
		//h stalokrokowe 0,12
		solution SDres1 = SD(x0, 0.12, epsilon, Nmax, limits);
		SDf_calls1 = solution::f_calls;
		SDg_calls1 = solution::g_calls;
		solution::clear_calls();

		cout << "CG 0,12" << endl;
		//h stalokrokowe 0,12
		solution CGres1 = CG(x0, 0.12, epsilon, Nmax, limits);
		CGf_calls1 = solution::f_calls;
		CGg_calls1 = solution::g_calls;
		solution::clear_calls();

		cout << "Newton 0,12" << endl;
		//h stalokrokowe 0,12
		solution Newtonres1 = Newton(x0, 0.12, epsilon, Nmax, limits);
		Newf_calls1 = solution::f_calls;
		Newg_calls1 = solution::g_calls;
		Newh_calls1 = solution::H_calls;
		solution::clear_calls();

		cout << "SD zm" << endl;
		//h zmiennokrokowe - bedzie obliczone zlotym podzialem
		solution SDres2 = SD(x0, -0.05, epsilon, Nmax, limits);
		SDf_calls2 = solution::f_calls;
		SDg_calls2 = solution::g_calls;
		solution::clear_calls();

		cout << "CG zm" << endl;
		//h zmiennokrokowe
		solution CGres2 = CG(x0, -0.05, epsilon, Nmax, limits);
		CGf_calls2 = solution::f_calls;
		CGg_calls2 = solution::g_calls;
		solution::clear_calls();

		cout << "Newton zm" << endl;
		//h zmiennokrokowe
		solution Newtonres2 = Newton(x0, -0.05, epsilon, Nmax, limits);
		Newf_calls2 = solution::f_calls;
		Newg_calls2 = solution::g_calls;
		Newh_calls2 = solution::H_calls;
		solution::clear_calls();
		*/
		/*f << x0(0)<<";"<<x0(1)<<";"<<SDres.x(0) << ";" << SDres.x(1) << ";" << SDres.y(0) << ";" << SDf_calls << ";" << SDg_calls << ";" << CGres.x(0)
				<< ";" << CGres.x(1) << ";" << CGres.y(0) << ";" << CGf_calls << ";" << CGg_calls << ";" <<
				Newtonres.x(0) << ";" << Newtonres.x(1) << ";" << Newtonres.y(0) << ";" << Newf_calls << ";" << Newg_calls
				<< ";" << Newh_calls << endl<<  ";"<<";" << SDres1.x(0) << ";" << SDres1.x(1) << ";" << SDres1.y(0) << ";" << SDf_calls1 << ";" << SDg_calls1 << ";"
				<< CGres1.x(0) << ";" << CGres1.x(1) << ";" << CGres1.y(0) << ";" << CGf_calls1 << ";" << CGg_calls1 << ";"
				<< Newtonres1.x(0) << ";" << Newtonres1.x(1) << ";" << Newtonres1.y(0) << ";" << Newf_calls1 << ";" <<
				Newg_calls1 << ";" << Newh_calls1 << endl  << ";"<<";" << SDres2.x(0) << ";" << SDres2.x(1) << ";" << SDres2.y(0) << ";" << SDf_calls2 << ";" << SDg_calls2 << ";"
				<< CGres2.x(0) << ";" << CGres2.x(1) << ";" << CGres2.y(0) << ";" << CGf_calls2 << ";" << CGg_calls2 << ";"
				<< Newtonres2.x(0) << ";" << Newtonres2.x(1) << ";" << Newtonres2.y(0) << ";" << Newf_calls2 << ";" << 
				Newg_calls2 << ";" << Newh_calls2 << endl;*/
		/*
		matrix x0(3, 1);
		x0(0) = x0(1) = x0(2) = 0.0;
		double krok = 0.01;
		double epsilon = 1 * 10 ^ -5;
		int Nmax = 20000;
		solution solCG = CG(x0, krok, epsilon, Nmax);
		cout << "Metoda CG z krokiem " << krok << endl;
		cout << "x0(0)= " << solCG.x(0) << endl;
		cout << "x0(1)= " << solCG.x(1) << endl;
		cout << "x0(2)= " << solCG.x(2) << endl;
		cout << "y= " << solCG.y(0) << endl;
		cout << "f_calls: " << solution::f_calls << endl;
		cout << "g_calls: " << solution::g_calls << endl;
		cout << "H_calls: " << solution::H_calls << endl;
		solution::clear_calls();
		*/
		solution X;
		matrix x0(3, 1);
		x0(0) = -297135;
		x0(1) = 700330;
		x0(2) = 0.00000107;
		X.x = x0;
		X.fit_fun();
#elif LAB_NO == 6
		static default_random_engine generate(unsigned(time(nullptr)));
		uniform_real_distribution<double> distribution(-10, 10);

		ofstream f("lab6_res.csv");
		ifstream start_point("x0.txt");
		double epsilon = 0.0001, w = 0.0;
		int Nmax = 1000;

		matrix limits(2, 3);
		limits(0, 0) = -10;
		limits(0, 1) = 10;
		limits(1, 0) = -10;
		limits(1, 1) = 10;
		limits(0, 2) = w;
		limits(1, 2) = 0;
		matrix x0(2, 1);

		for (int i = 0; i < 101; i++) {
			/*x0(0) = distribution(generate);
			x0(1) = distribution(generate);*/
			start_point >> x0(0);
			start_point >> x0(1);

			solution powell = Powell(x0, epsilon, Nmax, limits);
			f << x0(0) << ";" << x0(1) << ";" << powell.x(0) << ";" << powell.x(1) << ";" 
				<< powell.y(0) << ";" << powell.y(1) << ";" << solution::f_calls << endl;
			limits(0, 2) = limits(0, 2) + 0.01;
			solution::clear_calls();
		}
#endif
	}
	catch (char* EX_INFO)
	{
		std::cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
