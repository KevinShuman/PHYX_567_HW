// Created by Kevin Shuman on 1/24/20.
//  Homework2.c
//  
//
// The goal of this assignment is to numerically integrate a given function.
//
//
//
// Compile with...
// gcc -w -o Homework2 Homework2.c -lgsl -lgslcblas -lm
//
// Run with...
// ./Homework2


// Including the needed libraries
#include <stdio.h>
#include <math.h>

// This function takes in a values and computes the integrand. I named in Integrad on purpose.
double Integrad(double x);

int main()
{
    long long int i, j, N;
    
    // I allow YOU to give the number of steps you want.
    printf("Enter the step size:");
    scanf("%lli", &N);
    
    // Allocating memory for my arrays.
    double *x = (double *)malloc((size_t)sizeof(double)*(N+1));
    double *integrad = (double *)malloc((size_t)sizeof(double)*(N+1));
    
    // a is the values of the integral after each step we take. I suppose I could make this into an array
    // and have it flag at each order in step, but I've already have my results. Maybe an exercise for the reader.
    double a = 0.0;
    
    // Determining what values of x I want and plugging those into the Integrad function to create and array
    // I'll use in the integration process.
    for(i=1; i<=N; i++)
    {
        x[i] = -0.5 + ((double)(i)-0.5)*(1.0/((double)(N)));
        integrad[i] = Integrad(x[i]);
    }
    
    // I was curious as to what the integrand looked like, so this spits out a file with that info.
    /*
    FILE *out;
    out = fopen("Integrad.dat", "w");
    for(i=1; i<=N; i++) fprintf(out, "%.15e %.15e \n", x[i], integrad[i]);
    fclose(out);
     */
    
    // Setting up the f's for the integral.
    double result = 0.0;
    double f1, f2, f3, fn, fn1, fn2;
    f1 = 26.0/25.0*integrad[1];
    f2 = 21.0/24.0*integrad[2];
    f3 = 25.0/24.0*integrad[3];
    fn2 = 25.0/24.0*integrad[N-2];
    fn1 = 21.0/24.0*integrad[N-1];
    fn = 26.0/25.0*integrad[N];
    
    // Here's the summation. Notice I have i--. I'm actually starting from x = 0.5. My hopes were that
    // I might reduce some error later on by adding small things to small things at the beginning and
    // larger things to larger thing at the end, rather than the opposite by going the other way.
    // Didn't seem to do much, though.
    for(i=N; i>=1; i--)
    {
        if(i==1) a = f1;
        if(i==2) a = f2;
        if(i==3) a = f3;
        if(i==N-2) a = fn2;
        if(i==N-1) a = fn1;
        if(i==N) a = fn;
        if(i>=4 && i<=N-3)
        {
            a = integrad[i];
        }
        result += a/N;
    }
    
    // Gives the result.
    printf("\n\nThe result is: %.15e\n", result - 3.0);
    
    
    // Deallocating memory.
    free(x); free(integrad);
}


double Integrad(double x)
{
    double result, a, b, c, d;
    
    // Just as you'd expect for the integrand.
    a = 0.5 - x;
    b = pow(a,1.0/3.0);
    c = log(0.5 + x);
    //d = pow((0.5-x),-2.0/3.0);
    result = b/c+1.0/(b*b);
    
    return result;
}
