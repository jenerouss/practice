// Задание по практите Бехлера Никиты оптимизация метод деформируемого многогранника
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
using std::pow;
using std::sin;
using std::cos;
using std::fabs;
const double n = 2;
const double nr = n + 1;
double A = 1, B = 0.5, Y = 2, E = 0.2, K = 0;
double xl, xh, xt;
double xs = 0;
double o;
double ox = 0;
double xotr;
double f(double x)
{
	return pow(x, 2);
}
double fifthstep(double* x, double xt, double A, double xh)
{
	xotr = xt + A * (xt - xh);
	cout << "5func) xt= " << xt << " xh= " << xh << " xotr= " << xotr << " o= " << o << " xs= " << xs << " xl= " << xl << " ox= " << ox << endl;
	for (int i = 0; i < nr; i++)
	{
		cout << "5func)" << "x" << "[" << i << "] = " << x[i] << endl;
	}
	return xotr;
}
double fourthstep(double* x, double xt)
{
	ox = 0;
	for (int i = 0; i < nr; i++)
	{
		ox += f(x[i]) - f(xt);
	}
	o = pow((1 / (n + 1)) * pow(ox, 2), 0.5);
	if (o <= E)
	{
		cout << "o <= E, следовательно в качестве можно приближенного" <<
			"решения взять наилучшую точку текущего многогранника : " << xl << endl;
		return 0;
	}
	cout << "4func) xt= " << xt << " xh= " << xh << " xotr= " << xotr << " o= " << o << " xs= " << xs << " xl= " << xl << " ox= " << ox << endl;
	for (int i = 0; i < nr; i++)
	{
		cout << "4func)" << "x" << "[" << i << "] = " << x[i] << endl;
	}
	fifthstep(x, xt, A, xh);
	return o;
}
double thirdstep(double* x, double xh)
{
	xt = 0;
	for (int i = 0; i < nr; i++)
	{
		if (x[i] != xh)
			xt += x[i];
	}
	xt = (1 / n) * xt;
	cout << "3func) xt= " << xt << " xh= " << xh << " xotr= " << xotr << " o= " << o << " xs= " << xs << " xl= " << xl << " ox= " << ox << endl;
	for (int i = 0; i < nr; i++)
	{
		cout << "3func)" << "x" << "[" << i << "] = " << x[i] << endl;
	}
	fourthstep(x, xt);
	return xt;
}
double secondstep(double* x)
{
	for (int i = 0; i < nr; i++)
	{
		if (i >= 1)
		{
			if (min(x[i], x[i - 1]) < xl)
			{
				xl = min(x[i], x[i - 1]);
				//xlx = x[i];
				//xly = y[i];
			}
			if (max(x[i], x[i - 1]) > xh)
			{
				xh = max(x[i], x[i - 1]);
				//xhx = x[i];
				//xhy = y[i];
			}
		}
		else
		{
			xl = x[i];
			//xlx = x[i];
			//xly = y[i];
			xh = x[i];
			//xhx = x[i];
			//xhy = y[i];
		}
	}
	xs = 0;
	for (int i = 0; i < nr; i++)
	{
		if (x[i] > xl and x[i] < xh and x[i] > xs)
		{
			xs = x[i];
		}
	}
	cout << "2func) xt= " << xt << " xh= " << xh << " xotr= " << xotr << " o= " << o << " xs= " << xs << " xl= " << xl << " ox= " << ox << endl;
	for (int i = 0; i < nr; i++)
	{
		cout << "2func)" << "x" << "[" << i << "] = " << x[i] << endl;
	}
	thirdstep(x, xh);
	return xl, xh, xs;
}
int main()
{
	setlocale(0, "");
	srand(time(NULL));
	//cout << "Число Эпсилон(точность метода) Epsilon > 0 : ";
	//cin >> E;
	double xsz = 0;
	double xras = 0; // растяжение
	double* res = new double[nr];
	double* x = new double[nr];
	double* y = new double[nr];
	for (int i = 0; i < nr; i++)
	{
		x[i] = rand() % 50 + 1;
		y[i] = rand() % 50 + 1;
		while (x[i] == y[i])
		{
			y[i] = rand() % 50 + 1;
		}
		if (x[i] > y[i])
		{
			double h = x[i];
			x[i] = y[i];
			y[i] = h;
		}
		while (x[i] == x[i - 1] or x[i] == x[i - 2])
			x[i] = rand() % 50 + 1;
		res[i] = f(x[i]);
	}
	for (int i = 0; i < nr; i++)
	{
		cout << "x" << "[" << i << "] = " << x[i] << endl;
	}
	secondstep(x);
	while (o > E)
	//for (int i = 0; i < 10; i++)
	{
		cout << "o = " << o << endl;
		if (f(xotr) <= f(xl))
		{
			xras = xt + Y * (xotr - xt);
			cout << "1if" << endl;
			if (f(xras) < f(xl))
			{
				for (int j = 0; j < nr; j++)
				{
					if (x[j] == xh)
						x[j] = xras;
				}
				K = K + 1;
				secondstep(x);
			}
			else if (f(xras) >= f(xl))
			{
				for (int j = 0; j < nr; j++)
				{
					if (x[j] == xh)
						x[j] = xotr;
				}
				K = K + 1;
				secondstep(x);
			}
		}
		else if (f(xs) < f(xotr) and f(xotr) <= f(xh))
		{
			xsz = xt + B * (xh - xt);
			cout << "2elseif" << endl;
			for (int j = 0; j < nr; j++)
			{
				if (x[j] == xh)
					x[j] = xsz;
			}
			K = K + 1;
			secondstep(x);
		}
		else if (f(xl) < f(xotr) and f(xl) <= f(xs))
		{
			for (int j = 0; j < nr; j++)
			{
				if (x[j] == xh)
					x[j] = xotr;
			}
			cout << "3elseif" << endl;
			K = K + 1;
			secondstep(x);
		}
		else if (f(xotr) > f(xh))
		{
			for (int i = 0; i < nr; i++)
			{
				x[i] = xl + 0.5 * (x[i] - xl);
			}
			cout << "4elseif" << endl;
			K = K + 1;
			secondstep(x);
		}
	}
	cout << "r(xs)= " << f(xs) << endl << "r(xotr)= " << f(xotr) << endl << "r(xh)= " << f(xh) << endl << "r(xl)= " << f(xl) << endl;
	for (int i = 0; i < nr; i++)
	{
		cout << "x = " << x[i] << endl;
		//cout << "y = " << y[i] << endl;
		cout << "f(x" << i + 1 << ')' << '=' << res[i] << endl;
	}
	cout << "наилучшая: " << xl << endl;
	cout << "наихудшая: " << xh << endl;
	cout << "xs = " << xs << endl;
	cout << "центр тяжести: " << xt << endl;
	cout << "o = " << o << endl;
	cout << "отражение: " << f(xotr) << endl;
	cout << "растяжение: " << f(xras) << endl;
	cout << "сжатие: " << f(xsz) << endl;
	delete[] x;
	delete[] y;
	delete[] res;
	return 0;
}
