//
//  Homework3.c
//
// The aim of this assignment is to investigate stratified and non-stratified Monte Carlo integration for a shape whose boundary isn't
// one dimensional.
//
//  Created by Kevin Shuman on 2/10/20.
//

// Compile with
// gcc -w -o Homework3 Homework3.c -lgsl -lgslcblas -lm

// Including the libraries, where gsl_rng.h is that for random number generation
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int main()
{
    long int N, n, M, L, i, j, k, l, m;
    
    // I'm setting n to 1024. I'm using powers of 2 for both n and N since the stratification code needs evenly spaced boxes/
    n = pow(2,10);
    
    // L is the number of times we run the code to get our averages and standard deviation.
    L = 10;
    
    // M is related to the number of random points, N, used. I used 6 different values of N.
    M = 6;
    
    printf("The number of random points, N: %i \n", N);
    printf("The number of iterations,   n: %i \n", n);
    
    double **Zr, **Zi, *ZZ, *Std, *nStd;
    
    // I'm setting up the random number seed and defining variables for that here.
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    // More variables
    double scale, *Area, AreaB;
    FILE *out;
    
    // Allocating memory for arrays and zeroing them.
    Area = (double *)malloc((size_t)sizeof(double)*(L));
    nStd = (double *)malloc((size_t)sizeof(double)*(M));
    for(m=0; m<M; m++) nStd[m] = 0.0;
    Std = (double *)malloc((size_t)sizeof(double)*(M));
    for(m=0; m<M; m++) Std[m] = 0.0;
    
    double mean;
    
    // There are a lot of nested for loops here. The first is related to the number of points N.
    // This first set of loops is for the stratified Monte Carlo/
    for(m=0; m<M; m++)
    {
        // How I related N to M.
        N = pow(2,2*(m+3));
        
        // The next loop is to get different Monte Carlo runs for each N.
        for(l=0; l<L; l++)
        {
        
            // Allocating memory and zeroing arrays.
            Zr = (double **)malloc((size_t)sizeof(double)*(N));
            for(i=0; i<N; i++) Zr[i] = (double *)malloc((size_t)sizeof(double)*(n));
            Zi = (double **)malloc((size_t)sizeof(double)*(N));
            for(i=0; i<N; i++) Zi[i] = (double *)malloc((size_t)sizeof(double)*(n));
            ZZ = (double *)malloc((size_t)sizeof(double)*(N));
    
            for(i=0; i<N; i++)
            {
                for(j=0; j<n; j++)
                {
                    Zr[i][j] = 0.0;
                    Zi[i][j] = 0.0;
                }
            }

            
            // Here is how I stratify the total area. I have N points, so I need a box with sqrt(N) side lengths.
            // Each little box has a width called scale. Knowing this, I can randomly find a real and imaginary component
            // of my complex number in each box. These are the intital points, c talked about in the homework.
            k=0;
            scale = 4.0/(sqrt((double)(N)));
            for(i=0; i<sqrt(N); i++)
            {
                for(j=0; j<sqrt(N); j++)
                {
                    Zr[k][0] = scale*gsl_rng_uniform(r) - (2.0 - (double)(j)*scale);
                    Zi[k][0] = scale*gsl_rng_uniform(r) + (2.0 - (double)(i)*scale - scale);
                    if(Zr[k][0] >= 2.0 || Zi[k][0] >= 2.0 || Zr[k][0] <= -2.0 || Zi[k][0] <= -2.0) printf("These are out of bounds: %.6e %.6e %i %i %i \n", Zr[k][0], Zi[k][0], i, j, k);
                    k += 1;
                }
            }
    
            // I then use those c's here to find the count for each point, if there is one.
            // Once the count is reached, I call that value 20.0, since if I don't I get NaN's, which I can't plot.
            // I wasn't able to come up with a way to both stop at the count value and iterativly get my difference and Area,
            // which is why I have this problem.
            for(j=0; j<n-1; j++)
            {
                Area[l] = 0.0;
                for(i=0; i<N; i++)
                {
                    Zr[i][j+1] = Zr[i][j]*Zr[i][j] - Zi[i][j]*Zi[i][j] + Zr[i][0];
                    Zi[i][j+1] = 2.0*Zi[i][j]*Zr[i][j] + Zi[i][0];
                    ZZ[i] = sqrt(Zr[i][j]*Zr[i][j] + Zi[i][j]*Zi[i][j]);
                    if(ZZ[i] != ZZ[i]) ZZ[i] = 20.0;
                    if(ZZ[i]<2.0) Area[l] += 1.0;
                }
            }
    
            // Here I export a file for the Mandelbrot set. I use this to get my plot, but only for the largest values of M and L.
            if(m==M-1 && l==L-1)
            {
                out = fopen("Mandelbrot_C.txt","w");
                for(i=0; i<N; i++) fprintf(out, "%.6e %.6e %.6e\n", Zr[i][0], Zi[i][0], ZZ[i]);
                fclose(out);
            }
    
            // Freeing allocated memory
            for(i=0; i<N; i++)
            {
                free(Zr[i]); free(Zi[i]);
            }

            free(Zr); free(Zi); free(ZZ);
        }
    
        // Calculating the mean and standard deviation for each l in L.
        mean = 0.0;
        for(l=0; l<L; l++) mean += (Area[l]*16.0/N)/L;
        for(l=0; l<L; l++) Std[m] += sqrt(((Area[l]*16.0/N)-mean)*((Area[l]*16.0/N)-mean)/L);
    
        printf("The mean for N = %i points, stratified, is %f and the standard deviation is %f\n\n", N, mean, Std[m]);
    
        
        // Here I do the non-stratified Monte Carlo. The code is the same except for the random point generation and I gather my differences in area.
        for(l=0; l<L; l++)
        {
        
            Zr = (double **)malloc((size_t)sizeof(double)*(N));
            for(i=0; i<N; i++) Zr[i] = (double *)malloc((size_t)sizeof(double)*(n));
            Zi = (double **)malloc((size_t)sizeof(double)*(N));
            for(i=0; i<N; i++) Zi[i] = (double *)malloc((size_t)sizeof(double)*(n));
            ZZ = (double *)malloc((size_t)sizeof(double)*(N));
        
            for(i=0; i<N; i++)
            {
                for(j=0; j<n; j++)
                {
                    Zr[i][j] = 0.0;
                    Zi[i][j] = 0.0;
                }
            }
        
            // Here's where I generate my points, c.
            for(i=0; i<N; i++)
            {
                Zr[i][0] = 4.0*gsl_rng_uniform(r) - 2.0;
                Zi[i][0] = 4.0*gsl_rng_uniform(r) - 2.0;
                if(Zr[i][0] >= 2.0 || Zi[i][0] >= 2.0 || Zr[i][0] <= -2.0 || Zi[i][0] <= -2.0) printf("%.6e %.6e \n", Zr[i][0], Zi[i][0]);
            }
        
            AreaB = 0.0;
        
            // While I do the same things as before, I make sure to export the differences in area I get for each n.
            out = fopen("Area_n.dat", "w");
            for(j=0; j<n-1; j++)
            {
                Area[l] = 0.0;
                for(i=0; i<N; i++)
                {
                    Zr[i][j+1] = Zr[i][j]*Zr[i][j] - Zi[i][j]*Zi[i][j] + Zr[i][0];
                    Zi[i][j+1] = 2.0*Zi[i][j]*Zr[i][j] + Zi[i][0];
                    ZZ[i] = sqrt(Zr[i][j]*Zr[i][j] + Zi[i][j]*Zi[i][j]);
                    if(ZZ[i]<2.0) Area[l] += 1.0;
                }
                fprintf(out, "%i %.5e \n", j, (Area[l]-AreaB)*(16.0/N));
                AreaB = Area[l];
            }
            fclose(out);
        
            for(i=0; i<N; i++)
            {
                free(Zr[i]); free(Zi[i]);
            }
        
            free(Zr); free(Zi); free(ZZ);
        }
    
        mean = 0.0;
        for(l=0; l<L; l++) mean += (Area[l]*16.0/N)/L;
        for(l=0; l<L; l++) nStd[m] += sqrt(((Area[l]*16.0/N)-mean)*((Area[l]*16.0/N)-mean)/L);
    
        printf("The mean for N = %i points, not stratified, is %f and the standard deviation is %f\n\n", N, mean, nStd[m]);
    }
    
    // Once all is done, I grab my standard deviations for the stratified and non-stratified Monte Carlos and create a file to be plotted.
    out = fopen("Stand_Div.txt","w");
    for(m=0; m<M; m++) fprintf(out, "%i %f %f \n", (int)(pow(2,2*(m+3))), nStd[m], Std[m]);
    fclose(out);
    
    // Freeing more memory.
    free(Std); free(Area);
    gsl_rng_free (r);
}
