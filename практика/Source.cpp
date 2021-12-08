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
const double sx = 5;
const double sy = 2;
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
double* fifthstep(double** x,double* res,double* xl, double* xh, double* xs, double* xt,double ox, double o, double E, double A, double* xotr, double* xsz, double* xras)
{
	for (int i = 0; i < ny; i++)
	{
		xotr[i] = xt[i] + A * (xt[i] - xh[i]);
	}
	if (vivod)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				cout << "5func) x" << "[" << i << "]"
					<< " = (" << x[i][j] << ";" << x[i][j + 1] << ")" << endl;
				j++;
			}
			cout << "5func) res" << "[" << i << "]"
				<< " = " << res[i] << endl;
		}
		for (int i = 0; i < ny; i++)
		{
			cout << "5func) " << "xh = (" << xh[i] << ";" << xh[i + 1] << ") "
				<< "xs = (" << xs[i] << ";" << xs[i + 1] << ") "
				<< "xl = (" << xl[i] << ";" << xl[i + 1] << ") "
				<< "xt = (" << xt[i] << ";" << xt[i + 1] << ") "
				<< " ox = (" << ox << ") "
				<< " o = (" << o << ") "
				<< " xotr = (" << xotr[i] << ";" << xotr[i + 1] << ") "
				<< " xsz = (" << xsz[i] << ";" << xsz[i + 1] << ") "
				<< " xras = (" << xras[i] << ";" << xras[i + 1] << ")"
				<< endl;
			i++;
		}
	}
	return xotr;
}
double** fourthstep(double** x,double* res,double* xl, double* xh, double* xs, double* xt,double ox, double o, double E, double A, double* xotr, double* xsz, double* xras, double** saver)
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
			"решения взять наилучшую точку текущего многогранника : ";
		for (int i = 0; i < ny; i++)
		{
			cout << "(" << xl[i] << ";" << xl[i + 1] << ")";
			i++;
		}
		return 0;
	}
	else
	{
		cout << "o (" << o << ") > E (" << E << "), следовательно процесс поиска продолжается." << endl;
	}
	if (vivod)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				cout << "4func) x" << "[" << i << "]"
					<< " = (" << x[i][j] << ";" << x[i][j + 1] << ")" << endl;
				j++;
			}
			cout << "4func) res" << "[" << i << "]"
				<< " = " << res[i] << endl;
		}
		for (int i = 0; i < ny; i++)
		{
			cout << "4func) " << "xh = (" << xh[i] << ";" << xh[i + 1] << ") "
				<< "xs = (" << xs[i] << ";" << xs[i + 1] << ") "
				<< "xl = (" << xl[i] << ";" << xl[i + 1] << ") "
				<< "xt = (" << xt[i] << ";" << xt[i + 1] << ") "
				<< " ox = (" << ox << ") "
				<< " o = (" << o << ") "
				<< " xotr = (" << xotr[i] << ";" << xotr[i + 1] << ") "
				<< " xsz = (" << xsz[i] << ";" << xsz[i + 1] << ") "
				<< " xras = (" << xras[i] << ";" << xras[i + 1] << ")"
				<< endl;
			i++;
		}
	}
	for (int i = 4; i < sx; i++)
	{
		for (int j = 0; j < sy; j++)
		{
			saver[i][j] = o;
		}
	}
	fifthstep(x, res, xl, xh, xs, xt, ox, o, E, A, xotr, xsz, xras);
	return saver;
}
double** thirdstep(double** x, double* res, double* xl, double* xh, double* xs, double* xt, double ox, double o, double E, double A, double* xotr, double* xsz, double* xras, double** saver)
{
	double* sum = new double[ny] {};
	int k = 0;
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			if (x[i] != xh)
			{
				sum[k] += x[i][j];
			}
		}
		if (x[i] != xh)
		{
			k++;
			if (k >= ny)
				k--;
		}
	}
	for (int i = 0; i < ny; i++)
	{
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
				<< "xt = (" << xt[i] << ";" << xt[i + 1] << ") "
				<< " ox = (" << ox << ") "
				<< " o = (" << o << ") "
				<< " xotr = (" << xotr[i] << ";" << xotr[i + 1] << ") "
				<< " xsz = (" << xsz[i] << ";" << xsz[i + 1] << ") "
				<< " xras = (" << xras[i] << ";" << xras[i + 1] << ")"
				<< endl;
			i++;
		}
	}
	for (int i = 3; i < sx-1; i++)
	{
		for (int j = 0; j < sy; j++)
		{
			saver[i][j] = xt[j];
		}
	}
	fourthstep(x, res, xl, xh, xs, xt, ox, o, E, A, xotr, xsz, xras, saver);
	return saver;
}
double** secondstep(double** x, double* res, double* xl, double* xh, double* xs, double* xt, double ox, double o, double E, double A, double* xotr, double *xsz, double* xras, double** saver)
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			res[i] = f(x[i]);
		}
	}
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
				<< "xl = (" << xl[i] << ";" << xl[i + 1] << ") "
				<< "xt = (" << xt[i] << ";" << xt[i + 1] << ") "
				<< " ox = (" << ox << ") "
				<< " o = (" << o << ") "
				<< " xotr = (" << xotr[i] << ";" << xotr[i + 1] << ") "
				<< " xsz = (" << xsz[i] << ";" << xsz[i + 1] << ") "
				<< " xras = (" << xras[i] << ";" << xras[i + 1] << ")"
				<< endl;
			i++;
		}
	}
	for (int i = 0; i < sx-3; i++)
	{
		for (int j = 0; j < sy; j++)
		{
			saver[i][j] = xl[j];
		}
		i++;
		for (int j = 0; j < sy; j++)
		{
			saver[i][j] = xh[j];
		}
		i++;
		for (int j = 0; j < sy; j++)
		{
			saver[i][j] = xs[j];
		}
	}
	thirdstep(x, res, xl, xh, xs, xt, ox, o, E, A, xotr, xsz, xras, saver);
	return saver;
}
int main()
{
	setlocale(0, "");
	srand(time(NULL));
	//cout << "Число Эпсилон(точность метода) Epsilon > 0 : ";
	//cin >> E;
	double A = 1, B = 0.5, Y = 2, E = 0.2, K = 0;
	double ox = 0;
	double o = 0;
	double* res = new double[nx];
	double* xl = new double[ny];
	double* xh = new double[ny];
	double* xt = new double[ny] {};
	double* xs = new double[ny];
	double* xotr = new double[ny] {};
	double* xras = new double[ny] {};
	double* xsz = new double[ny] {};
	double** saver;
	saver = new double* [sx];
	double** x;
	x = new double* [nx];
	for (int i = 0; i < sx; i++)
	{
		saver[i] = new double[sy];
	}
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
	saver = secondstep(x,res,xl,xh,xs,xt,o,ox,E,A,xotr,xsz,xras,saver);
	for (int i = 0; i < sx; i++)
	{
		for (int j = 0; j < sy; j++)
		{
			xl[j] = saver[i][j];
		}
		i++;
		for (int j = 0; j < sy; j++)
		{
			xh[j] = saver[i][j];
		}
		i++;
		for (int j = 0; j < sy; j++)
		{
			xs[j] = saver[i][j];
		}
		i++;
		for (int j = 0; j < sy; j++)
		{
			xt[j] = saver[i][j];
		}
		i++;
		for (int j = 0; j < sy; j++)
		{
			o = saver[i][j];
		}
	}
	//while (o > E)
	for (int j = 0; j < 10; j++)
	{
		for (int i = 0; i < sx; i++)
		{
			for (int j = 0; j < sy; j++)
			{
				xl[j] = saver[i][j];
			}
			i++;
			for (int j = 0; j < sy; j++)
			{
				xh[j] = saver[i][j];
			}
			i++;
			for (int j = 0; j < sy; j++)
			{
				xs[j] = saver[i][j];
			}
			i++;
			for (int j = 0; j < sy; j++)
			{
				xt[j] = saver[i][j];
			}
			i++;
			for (int j = 0; j < sy-1; j++)
			{
				o = saver[i][j];
			}
		}
		if (vivod)
			cout << "o = " << o << endl;
		if (f(xotr) <= f(xl))
		{
			if (vivod)
				cout << "1if" << endl;
			for (int i = 0; i < ny; i++)
			{
				xras[i] = xt[i] + Y * (xotr[i] - xt[i]);
			}
			if (f(xras) < f(xl))
			{
				if (vivod)
					cout << "1ifif" << endl;
				for (int i = 0; i < nx; i++)
				{
					if (*x[i] == *xh)
					{
						x[i] = xras;
					}
				}
				K = K + 1;
				secondstep(x, res, xl, xh, xs, xt, o, ox, E, A, xotr, xsz, xras,saver);
			}
			else if (f(xras) >= f(xl))
			{
				if (vivod)
					cout << "1elseif" << endl;
				for (int i = 0; i < nx; i++)
				{
					if (*x[i] == *xh)
					{
						x[i] = xotr;
					}
				}
				K = K + 1;
				secondstep(x, res, xl, xh, xs, xt, o, ox, E, A, xotr, xsz, xras,saver);
			}
		}
		else if (f(xs) < f(xotr) and f(xotr) <= f(xh))
		{
			if (vivod)
				cout << "2elseif" << endl;
			for (int i = 0; i < ny; i++)
			{
				xsz[i] = xt[i] + B * (xh[i] - xt[i]);
			}
			for (int i = 0; i < nx; i++)
			{
				if (*x[i] == *xh)
					x[i] = xsz;
			}
			K = K + 1;
			secondstep(x, res, xl, xh, xs, xt, o, ox, E, A, xotr, xsz, xras,saver);
		}
		else if (f(xl) < f(xotr) and f(xl) <= f(xs))
		{
			if (vivod)
				cout << "3elseif" << endl;
			for (int i = 0; i < nx; i++)
			{
				if (*x[i] == *xh)
					x[i] = xotr;
			}
			K = K + 1;
			secondstep(x, res, xl, xh, xs, xt, o, ox, E, A, xotr, xsz, xras,saver);
		}
		else if (f(xotr) > f(xh))
		{
			if (vivod)
				cout << "4elseif" << endl;
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					x[i][j] = xl[j] + 0.5 * (x[i][j] - xl[j]);
				}
			}
			K = K + 1;
			secondstep(x, res, xl, xh, xs, xt, o, ox, E, A,xotr,xsz, xras,saver);
		}
	}
	delete[] res;
	delete[] xl;
	delete[] xh;
	delete[] xt;
	delete[] xs;
	delete[] xotr;
	delete[] xras;
	delete[] xsz;
	for (int i = 0; i < nx; i++)
	{
		delete[] x[i];
	}
	delete[] x;
	for (int i = 0; i < sx; i++)
	{
		delete[] saver[i];
	}
	delete[] saver;
	return 0;
}
