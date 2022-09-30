// Compile with
// gcc -o Homework_5 Homework_5.c -lm -lgslcblas -lgsl
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

// I have two functions, one that is the differential equation and the other is the Jacobian, which I calculate by hand.
int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

int main (void)
{	
	int i, j, k;
	
	// N is the number of points for the E0 and omega parameters, so N*N is the total resolution of chaotic region diagram.
	// Here I choose N to be 20, but I have differernt data sets I have calculated. Two with N = 100 and one with N = 500.
	int N = 100;
	
	// I'm allocating memory here.
	double *E0_array = (double *)malloc((size_t)sizeof(double) * (N));
	double *omega_array = (double *)malloc((size_t)sizeof(double) * (N));
	int *sensitivity = (int *)malloc((size_t)sizeof(int) * (N*N));
	
	// The way I calculate the histogram is with this sensitivity array. I sent all the values to zero, but if the tajectories
	// seperate a given distance, I call the tolerance, it switchs the value to one.
	for(i=0; i<N*N; i++) sensitivity[i] = 0;
	
	// Here is the tolerance I mentioned above.
	double tolerance = 0.01;
	
	// Here I'm creating my grid in parameter space I will go over. It is even spaced in both directions.
	for(i=0; i<N; i++)
	{
		E0_array[i] = (double)(i)*1.0e6/N;
		omega_array[i] = (double)(i)*30.0/N;
	}

	// Defining my variables.
	double g;
	double mu;
	double Q;
	double E0;
	double M;
	double L;
	double omega;
	
	double *params1;
	double *params2;
	
	double T1 = 0.0; double t11 = 30.0;
	double T2 = 0.0; double t12 = 30.0;
	double y1[2] = { 0.0, 0.0 };
	double y2[2] = { 0.0001, 0.0 };
	
	FILE *out;
	
	double ti1;
	double ti2;
	int status1;
	int status2;
	
	
	// This part I calculate the tajectories for the case Charles asked us to do; it's the same code,
	// but with out the sensitivity stuff.
	g = 9.81;
	mu = 0.1;
	Q = 0.00000001;
	E0 = 1000000.0;
	M = 0.0005;
	L = 0.1;
	omega = 11.0;
	
	params1 = (double *)malloc((size_t)sizeof(double) * (7));
	params2 = (double *)malloc((size_t)sizeof(double) * (7));
	
	params1[0] = g;
	params1[1] = mu;
	params1[2] = Q;
	params1[3] = E0;
	params1[4] = M;
	params1[5] = L;
	params1[6] = omega;
	
	params2[0] = g;
	params2[1] = mu;
	params2[2] = Q;
	params2[3] = E0;
	params2[4] = M;
	params2[5] = L;
	params2[6] = omega;
	
	gsl_odeiv2_system sys1 = {func, jac, 2, params1};
	gsl_odeiv2_driver *d1 = gsl_odeiv2_driver_alloc_y_new (&sys1, gsl_odeiv2_step_rk4imp, 1e-6, 1e-6, 0.0);
	
	gsl_odeiv2_system sys2 = {func, jac, 2, params2};
	gsl_odeiv2_driver *d2 = gsl_odeiv2_driver_alloc_y_new (&sys2, gsl_odeiv2_step_rk4imp, 1e-6, 1e-6, 0.0);


	T1 = 0.0; t11 = 30.0;
	T2 = 0.0; t12 = 30.0;
			
	y1[0] = 0.0, y1[1] = 0.0;
	y2[0] = 0.0001; y2[1] = 0.0;
			
	out = fopen("Trajectories.dat","w");
	for (k = 1; k <= 300; k++)
	{
		ti1 = (double)(k) * t11 / 300.0;
		ti2 = (double)(k) * t12 / 300.0;
		status1 = gsl_odeiv2_driver_apply (d1, &T1, ti1, y1);
		status2 = gsl_odeiv2_driver_apply (d2, &T2, ti2, y2);
					
		if (status1 != GSL_SUCCESS || status2 != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status1);
			printf ("error, return value=%d\n", status2);
			break;
		}
		// I create a dat file of the trajectories.
		fprintf(out, "%.6e %.6e %.6e \n", ti1, y1[0], y2[0]);
	}
	fclose(out);

	free(params1); free(params2);
	gsl_odeiv2_driver_free (d1); gsl_odeiv2_driver_free (d2);
	
	
	
	// Using my grid of parameters, I run through the list point by point and solving the ODE's for the two initial condiditions.
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			g = 9.81;
			mu = 0.1;
			Q = 0.00000001;
			E0 = E0_array[i];
			M = 0.0005;
			L = 0.1;
			omega = omega_array[j];
	
			// The way the ODE solver works is that I have to solve the ODEs with different initial conditions completely seperately.
			params1 = (double *)malloc((size_t)sizeof(double) * (7));
			params2 = (double *)malloc((size_t)sizeof(double) * (7));
	
			params1[0] = g;
			params1[1] = mu;
			params1[2] = Q;
			params1[3] = E0;
			params1[4] = M;
			params1[5] = L;
			params1[6] = omega;
	
			params2[0] = g;
			params2[1] = mu;
			params2[2] = Q;
			params2[3] = E0;
			params2[4] = M;
			params2[5] = L;
			params2[6] = omega;
	
			// Here I'm telling my ODE solver what ODE I'm solving, giving it my Jacobian, and handing it my parameters.
			gsl_odeiv2_system sys1 = {func, jac, 2, params1};
			
			// This is a wrapper that allocates all the memory I need and creates my workspace. "gsl_odeiv2_step_rk4imp"
			// in the function argument is the type of ODE solve I am using. In this case I'm using an implicit Gaussian 4th order Runge-Kutta.
			gsl_odeiv2_driver *d1 = gsl_odeiv2_driver_alloc_y_new (&sys1, gsl_odeiv2_step_rk4imp, 1e-6, 1e-6, 0.0);
	
			// I do the same for the other initial condition.
			gsl_odeiv2_system sys2 = {func, jac, 2, params2};
			gsl_odeiv2_driver *d2 = gsl_odeiv2_driver_alloc_y_new (&sys2, gsl_odeiv2_step_rk4imp, 1e-6, 1e-6, 0.0);


			// Here I'm giving the initial and final times.
			T1 = 0.0; t11 = 30.0;
			T2 = 0.0; t12 = 30.0;
			
			// Here I'm giving the initial conditions.
			y1[0] = 0.0, y1[1] = 0.0;
			y2[0] = 0.0001; y2[1] = 0.0;

			// This is where I solve the ODE. I create the time I'm stepping to and pass the ODE solver the time I am at and the new time.
			// It then changes time where I'm at to the new time and repeats.
			
			// Below is were I determine if my tajectories are seperating. If they seperate far enough apart, we break the for loop
			// and go onto the next point in parameter space.
			for (k = 1; k <= 300; k++)
			{
				ti1 = (double)(k) * t11 / 300.0;
				ti2 = (double)(k) * t12 / 300.0;
				status1 = gsl_odeiv2_driver_apply (d1, &T1, ti1, y1);
				status2 = gsl_odeiv2_driver_apply (d2, &T2, ti2, y2);
					
				if (status1 != GSL_SUCCESS || status2 != GSL_SUCCESS)
				{
					printf ("error, return value=%d\n", status1);
					printf ("error, return value=%d\n", status2);
					break;
				}
				
				if(abs(y1[0]-y2[0])>tolerance)
				{
					sensitivity[i*N + j] = 1;
					break;
				}
			}
			// It's just nice to see where the code is in its calculation.
			printf("%i %i\n", i, j);
			
			// Freeing allocated memory.
			free(params1); free(params2);
			gsl_odeiv2_driver_free (d1); gsl_odeiv2_driver_free (d2);
		}
	}
	
	
	// I create a file with al the info I need to make a histogram of the data.
	// To make this histogram I use Mathematica (I would use Matlab, but I don't have it on my laptop).
	out = fopen("Hist_Data_00_1_short.dat","w");
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++) fprintf(out, "%.6e %.6e %i \n", E0_array[i], omega_array[j], sensitivity[i*N+j]);
	}
	fclose(out);
	
	return 0;
}

// This is the function that describes the ODE.
int func (double t, const double y[], double f[], void *params)
{
	double *param = (double *)params;
	double g = param[0];
	double mu = param[1];
	double Q = param[2];
	double E0 = param[3];
	double M = param[4];
	double L = param[5];
	double omega = param[6];
	
	// Because it is second order, I need to seperate it into two first order equations.
	// The ODE solver is picky in that is needs yoo to have the first f as the trivial ODE.
	f[0] = y[1];
	f[1] = (g/2.0)*sin(y[0]) - mu*y[1] + (Q*E0/(M*L))*cos(y[0])*cos(omega*t);
	
	free(param);
	return GSL_SUCCESS;
}

// Here's my Jacobian function.
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	double *param = (double *)params;
	double g = param[0];
	double mu = param[1];
	double Q = param[2];
	double E0 = param[3];
	double M = param[4];
	double L = param[5];
	double omega = param[6];
	
	// Somehow this works. I took this part of the code for the example code on the website https://www.gnu.org/software/gsl/doc/html/ode-initval.html
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix;
	
	// Here are the terms of my Jacobian.
	gsl_matrix_set (m, 1, 1, (g/2.0)*cos(y[0])-(Q*E0/(M*L))*sin(y[0])*cos(omega*t));
	gsl_matrix_set (m, 1, 0, -mu);
	gsl_matrix_set (m, 0, 1, 0.0);
	gsl_matrix_set (m, 0, 0, 1.0);
	
	// The ODE solver also requires me to take any time derivatives of the f's
	dfdt[0] = -omega*(Q*E0/(M*L))*cos(y[0])*sin(omega*t);
	dfdt[1] = 0.0;
	
	free(param);
	return GSL_SUCCESS;
}