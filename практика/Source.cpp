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
const double nx = n + 1;
const double ny = 2;
bool vivod = true;
double f(double* x)
{
	double* d = new double[ny];
	for (int j = 0; j < nx; j++)
	{
		d[j] = x[j];
	}
	return (pow(d[0], 2) + pow(d[1],2));
}
/*double fifthstep(double** x,double* res,double* xl, double* xh, double* xs, double* xt)
{
	xotr = xt + A * (xt - xh);
	if (vivod)
	{
		cout << "5func) xt= " << xt << " xh= " << xh << " xotr= " << xotr << " o= " << o << " xs= " << xs << " xl= " << xl << " ox= " << ox << endl;
		for (int i = 0; i < nx; i++)
		{
			cout << "5func)" << "x" << "[" << i << "] = " << x[i] << endl;
		}
	}
	return xotr;
}
double fourthstep(double** x,double* res,double* xl, double* xh, double* xs, double* xt)
{
	ox = 0;
	for (int i = 0; i < nx; i++)
	{
		ox += f(x[i]) - f(xt);
	}
	o = pow((1 / (n + 1)) * pow(ox, 2), 0.5);
	if (o <= E)
	{
		cout << "o (" << o << ") <= E (" << E << "), следовательно в качестве можно приближенного" <<
			"решения взять наилучшую точку текущего многогранника : " << xl << endl;
		return 0;
	}
	if (vivod)
	{
		cout << "4func) xt= " << xt << " xh= " << xh << " xotr= " << xotr << " o= " << o << " xs= " << xs << " xl= " << xl << " ox= " << ox << endl;
		for (int i = 0; i < nx; i++)
		{
			cout << "4func)" << "x" << "[" << i << "] = " << x[i] << endl;
		}
	}
	fifthstep(x, xt, A, xh);
	return o;
}*/
double* thirdstep(double** x, double* res, double* xl, double* xh, double* xs, double* xt)
{
	double* sum = new double[ny] {};
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			if (x[i] != xh)
			{
				cout << "if " << x[i][j] << endl;
				sum[i] += x[i][j];
				cout << sum[i] << endl;
			}
		}
		xt[i] = (1 / n) * sum[i];
	}
	delete[] sum;
	if (vivod)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				cout << "3func) x" << "[" << i << "]"
					<< " = (" << x[i][j] << ";" << x[i][j + 1] << ")" << endl;
				j++;
			}
			cout << "3func) res" << "[" << i << "]"
				<< " = " << res[i] << endl;
		}
		for (int i = 0; i < ny; i++)
		{
			cout << "3func) " << "xh = (" << xh[i] << ";" << xh[i + 1] << ") "
				<< "xs = (" << xs[i] << ";" << xs[i + 1] << ") "
				<< "xl = (" << xl[i] << ";" << xl[i + 1] << ") "
				<< "xt = (" << xt[i] << ";" << xt[i + 1] << ")"
				<< endl;
			i++;
		}
	}
	//fourthstep(x, xt);
	return xt;
}
double* secondstep(double** x,double* res,double* xl, double* xh, double* xs, double* xt)
{
	for (int i = 0; i < nx; i++)
	{
		if (i >= 1)
		{
			if (min(res[i], res[i - 1]) < f(xl))
			{
				xl = x[i];
			}
			if (max(res[i], res[i - 1]) > f(xh))
			{
				xh = x[i];
			}
		}
		else
		{
			xl = x[i];
			xh = x[i];
		}
		xs = x[i];
	}
	for (int i = 0; i < nx; i++)
	{
		if (res[i] > f(xl) and res[i] < f(xh) and res[i] > xs[i])
		{
			xs = x[i];
		}
	}
	if (vivod)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				cout << "2func) x" << "[" << i << "]"
					<< " = (" << x[i][j] << ";" << x[i][j + 1] << ")" << endl;
				j++;
			}
			cout << "2func) res" << "[" << i << "]"
				<< " = " << res[i] << endl;
		}
		for (int i = 0; i < ny; i++)
		{
			cout << "2func) " << "xh = (" << xh[i] << ";" << xh[i+1] << ") " 
				<< "xs = (" << xs[i] << ";" << xs[i + 1] << ") "
				<< "xl = (" << xl[i] << ";" << xl[i + 1] << ")" << endl;
			i++;
		}
	}
	thirdstep(x, res, xl, xh, xs, xt);
	return xl, xh, xs;
}
int main()
{
	setlocale(0, "");
	srand(time(NULL));
	//cout << "Число Эпсилон(точность метода) Epsilon > 0 : ";
	//cin >> E;
	double xsz = 0;
	double xras = 0;
	double A = 1, B = 0.5, Y = 2, E = 0.2, K = 0;
	double o;
	double ox = 0;
	double* res = new double[nx];
	double* xl = new double[ny];
	double* xh = new double[ny];
	double* xt = new double[ny];
	double* xs = new double[ny];
	double* xotr = new double[ny];
	double** x;
	x = new double* [nx];
	for (int i = 0; i < nx; i++) 
	{
		x[i] = new double[ny];
	}
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			x[i][j] = rand() % 50 + 1;
			while (x[i][j] == x[i][j - 1])
			{
				x[i][j] = rand() % 50 + 1;
			}
			res[i] = f(x[i]);
		}
	}
	if (vivod)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				cout << "x" << "[" << i << "]" << " = (" << x[i][j] << ";" << x[i][j+1] << ")" << endl;
				j++;
			}
			cout << "res" << "[" << i << "]" << " = " << res[i] << endl;
		}
	}
	secondstep(x,res,xl,xh,xs,xt);
	/*while (o > E)
	{
		if (vivod)
			cout << "o = " << o << endl;
		if (f(xotr) <= f(xl))
		{
			xras = xt + Y * (xotr - xt);
			if (vivod)
				cout << "1if" << endl;
			if (f(xras) < f(xl))
			{
				for (int j = 0; j < nx; j++)
				{
					if (x[j] == xh)
						x[j] = xras;
				}
				K = K + 1;
				secondstep(x);
			}
			else if (f(xras) >= f(xl))
			{
				for (int j = 0; j < nx; j++)
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
			if (vivod)
				cout << "2elseif" << endl;
			for (int j = 0; j < nx; j++)
			{
				if (x[j] == xh)
					x[j] = xsz;
			}
			K = K + 1;
			secondstep(x);
		}
		else if (f(xl) < f(xotr) and f(xl) <= f(xs))
		{
			for (int j = 0; j < nx; j++)
			{
				if (x[j] == xh)
					x[j] = xotr;
			}
			if (vivod)
				cout << "3elseif" << endl;
			K = K + 1;
			secondstep(x);
		}
		else if (f(xotr) > f(xh))
		{
			for (int i = 0; i < nx; i++)
			{
				x[i] = xl + 0.5 * (x[i] - xl);
			}
			if (vivod)
				cout << "4elseif" << endl;
			K = K + 1;
			secondstep(x);
		}
	}
	if (vivod)
	{
		cout << "r(xs)= " << f(xs) << endl << "r(xotr)= " << f(xotr) << endl << "r(xh)= " << f(xh) << endl << "r(xl)= " << f(xl) << endl;
	for (int i = 0; i < nx; i++)
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
	}*/
	for (int i = 0; i < nx; i++)
	{
		delete[] x[i];
	}
	delete[] x;
	delete[] res;
	delete[] xl;
	delete[] xh;
	delete[] xt;
	delete[] xs;
	delete[] xotr;
	return 0;
}
