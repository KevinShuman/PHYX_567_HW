// Compile with
// gcc -o Homework_6 Homework_6.c -lm -lgslcblas -lgsl
//
// The goal of this homework is to investigate the sensitivity to initial conditions for a given differential equation.

// The libraries I am using.
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// These are the libraries for the ODE solver.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define PI 3.14159265358979323846

// I have two functions, one that is the differential equation and the other is the Jacobian, which I calculate by hand.
int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

int main (void)
{	
	int i, j, k;
	
	// N is the number of time steps I take, not including those from the adaptive steping.
	int N = 1000;

	// Defining my variables.
	double A;
	double ro;
	double omega;
	double m;
	double eta;
	double D;
	double a;
	
	double *params;
	
	double T; 
	double t;
	double y[4] = {0.0, 1500.0, 0.0, 0.0};
	
	FILE *out;
	
	double ti;
	int status;
	
	// I calculate the tajectories.
	params = (double *)malloc((size_t)sizeof(double) * (7));
	
	// My attempt at varying the diameter.
	D = 5.94e-6;
	//D = 10.0e-6;
	//D = 100.0e-6;
	//D = 1000.0e-6;
	//D = 10000.0e-6;
	
	A = PI*D*D*0.25;
	m = 4000.0*PI*D*D*D/24.0;
	ro = 1.0;
	omega = 10.0;
	eta = 2.0e-5;
	a = 1000.0;
	
	params[0] = A;
	params[1] = ro;
	params[2] = omega;
	params[3] = m;
	params[4] = eta;
	params[5] = D;
	params[6] = a;
	
	// Setting up the ODE solver. I use a "variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form"
	// here, but I have attempted all the other solvers without any success.
	gsl_odeiv2_system sys = {func, jac, 4, params};
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-4, 1e-4, 0.0);
	
	// Define my initial time and my final time.
	T = 0.0;
	t = 2.0*PI*a/omega;
	
	y[0] = 0.0; y[1] = 1500.0; y[2] = 0.0; y[3] = 0.0;
	
	// I run through all my steps with the stepperand save the tajectory to file.
	out = fopen("Trajectories_594.dat","w");
	for (k = 0; k <= N; k++)
	{
		ti = (double)(k) * t / ((double)(N));
		
		// Uncomment to have it tell you what step you are on while calculating, if that's your business.
		printf("Step: %i\n", k);
		status = gsl_odeiv2_driver_apply (d, &T, ti, y);
		
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		if (y[1]<= 0.0) break;
		
		// I create a dat file of the trajectory and velocities at each time.
		fprintf(out, "%.6e %.6e %.6e %.6e %.6e \n", ti, y[0], y[1], y[2], y[3]);
	}
	fclose(out);

	
	free(params);
	gsl_odeiv2_driver_free (d);

	return 0;
}

// This is the function that describes the ODE.
int func (double t, const double y[], double f[], void *params)
{
	double x, Y, u, v, wx, wy, V, Cx, Sx, Cy, Sy;

	double *param = (double *)params;
	double ro = param[1];
	double omega = param[2];
	double m = param[3];
	double eta = param[4];
	double d = param[5];
	double a = param[6];
	
	double A = PI*d*d*0.25;
	
	//(void)(t);
	x = y[0];
	Y = y[1];
	u = y[2];
	v = y[3];
	
	Cx = cos(PI*x/(2.0*a));
	Sx = sin(PI*x/(2.0*a));
	Cy = cos(PI*Y/(2.0*a));
	Sy = sin(PI*Y/(2.0*a));
	
	wx = -(PI*omega/2.0)*Cx*Cy;
	wy = -(PI*omega/2.0)*Sx*Sy;
	
	V = sqrt((wx-u)*(wx-u) + (wy-v)*(wy-v));
	
	// Because it is second order, I need to seperate it into two first order equations.
	// The ODE solver is picky in that is needs you to have the first f as the trivial ODE.
	
	f[0] = y[2];
	f[1] = y[3];
	f[2] =         (12.0*A*eta/(m*d))*(1.0 - (1.0/6.0)*pow(ro*V*d/eta, 1.5))*(wx - u);
	f[3] = -10.0 + (12.0*A*eta/(m*d))*(1.0 - (1.0/6.0)*pow(ro*V*d/eta, 1.5))*(wy - v);
	
	//printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e \n", t, f[0], f[1], f[2], f[3], y[0], y[1], y[2], y[3]);
	
	free(param);
	return GSL_SUCCESS;
}

// Here's my Jacobian function.
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	double x, Y, u, v, wx, wy, V, Cx, Sx, Cy, Sy;
	double duwxu, dvwxu, dxwxu, dywxu, duwyv, dvwyv, dxwyv, dywyv;
	double duV, dvV, dxV, dyV;
	double Jxu, Jxv, Jxx, Jxy;
	double Jyu, Jyv, Jyx, Jyy;
	
	double *param = (double *)params;
	double ro = param[1];
	double omega = param[2];
	double M = param[3];
	double eta = param[4];
	double d = param[5];
	double a = param[6];
	
	double A = PI*d*d*0.25;
	
	//(void)(t);
	u = y[0];
	v = y[1];
	x = y[2];
	Y = y[3];
	
	Cx = cos(PI*x/(2.0*a));
	Sx = sin(PI*x/(2.0*a));
	Cy = cos(PI*Y/(2.0*a));
	Sy = sin(PI*Y/(2.0*a));
	
	wx = -(PI*omega/2.0)*Cx*Cy;
	wy = -(PI*omega/2.0)*Sx*Sy;
	
	V = sqrt((wx-u)*(wx-u) + (wy-v)*(wy-v));
	
	duwxu = -1.0;
	dvwxu = 0.0;
	dxwxu = (omega*PI*PI/(4.0*a))*Sx*Cy;
	dywxu = (omega*PI*PI/(4.0*a))*Cx*Sy;
	
	duwyv = 0.0;
	dvwyv = -1.0;
	dxwyv = -(omega*PI*PI/(4.0*a))*Cx*Sy;
	dywyv = -(omega*PI*PI/(4.0*a))*Sx*Cy;
	
	duV = (u-wx)/V;
	dvV = (v-wy)/V;
	dxV = (PI*PI*omega/(4.0*a*V))*(Sx*Cy*(wx-u) - Cx*Sy*(wy-v));
	dyV = (PI*PI*omega/(4.0*a*V))*(Cx*Sy*(wx-u) - Sx*Cy*(wy-v));
	
	
	double arg1 = 24.0*eta/(ro*d);
	double arg2 = 4.0*sqrt(ro*d*V/eta)*V;
	double arg3 = 6.0*sqrt(ro*d*V/eta);
	
	Jxu = (A*ro/(2.0*M))*( ( arg1 - arg2 )*duwxu -arg3*(wx-u)*duV );
	Jxv = (A*ro/(2.0*M))*( ( arg1 - arg2 )*dvwxu -arg3*(wx-u)*dvV );
	Jxx = (A*ro/(2.0*M))*( ( arg1 - arg2 )*dxwxu -arg3*(wx-u)*dxV );
	Jxy = (A*ro/(2.0*M))*( ( arg1 - arg2 )*dywxu -arg3*(wx-u)*dyV );
	
	Jyu = (A*ro/(2.0*M))*( ( arg1 - arg2 )*duwyv -arg3*(wy-v)*duV );
	Jyv = (A*ro/(2.0*M))*( ( arg1 - arg2 )*dvwyv -arg3*(wy-v)*dvV );
	Jyx = (A*ro/(2.0*M))*( ( arg1 - arg2 )*dxwyv -arg3*(wy-v)*dxV );
	Jyy = (A*ro/(2.0*M))*( ( arg1 - arg2 )*dywyv -arg3*(wy-v)*dyV );
	
	// Somehow this works. I took this part of the code for the example code on the website https://www.gnu.org/software/gsl/doc/html/ode-initval.html
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
	gsl_matrix * m = &dfdy_mat.matrix;
	
	//printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e \n", t, Jxx, Jxy, Jxu, Jxv, Jyx, Jyy, Jyu, Jyv);
	
	// Here are the terms of my Jacobian.
	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 0.0);
	gsl_matrix_set (m, 0, 2, 1.0);
	gsl_matrix_set (m, 0, 3, 0.0);
	gsl_matrix_set (m, 1, 0, 0.0);
	gsl_matrix_set (m, 1, 1, 0.0);
	gsl_matrix_set (m, 1, 2, 0.0);
	gsl_matrix_set (m, 1, 3, 1.0);
	gsl_matrix_set (m, 2, 0, Jxx);
	gsl_matrix_set (m, 2, 1, Jxy);
	gsl_matrix_set (m, 2, 2, Jxu);
	gsl_matrix_set (m, 2, 3, Jxv);
	gsl_matrix_set (m, 3, 0, Jyx);
	gsl_matrix_set (m, 3, 1, Jyy);
	gsl_matrix_set (m, 3, 2, Jyu);
	gsl_matrix_set (m, 3, 3, Jyv);
	
	
	// The ODE solver also requires me to take any time derivatives of the f's
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;
	
	free(param);
	return GSL_SUCCESS;
}