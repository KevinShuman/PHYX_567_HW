//
//  04_kevin_shuman.c
//  
//
// The goal of this code is to take a WAV file with some noisy music, and reduce the
// noise using and optimal filter.
//
//
//  Created by Kevin Shuman on 2/21/20.
//

// Compile with
// gcc -o 04_kevin_shuman 04_kevin_shuman.c -lgsl -lgslcblas -lm -lsndfile

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_fft_real.h> // GSL libraries for fft's
#include <gsl/gsl_fft_halfcomplex.h>
#include <string.h>

// The following header file comes from a library:
// http://www.mega-nerd.com/libsndfile/
// I used Homebrew to install it on my Mac using the following command
// brew install libsndfile
// Some permission issues came up, but following what the errors tell you to do
// solves those problems
#include <sndfile.h>

// This is a function I stole from https://stackoverflow.com/questions/8424170/1d-linear-convolution-in-ansi-c-code
void convolve(double *Signal, int SignalLen, double *Kernel, int KernelLen, double *Result);

int main()
{
    // The first step is to read in the song and noise
    // I use the libsndfile library to do this.
    // There is a code Input_WAVE_Test.c that gives an
    // example as to how to do this.
    
    SNDFILE *tango, *noise; // These are specific structures for the library
    SF_INFO info_t, info_n;
    
    int num_t_channels, num_n_channels;
    int num_t, num_n, num_t_items, num_n_items;
    int *buf_t_i, *buf_n_i;
    double *buf_t_1, *buf_t_2;
    double *buf_t, *buf_n;
    int f_t, f_n, sr_t, sr_n, c_t, c_n;
    int i, j;
    
    // Open the WAV files
    info_t.format = 0;
    info_n.format = 0;
    tango = sf_open("Lacumparsita-noisy.wav",SFM_READ,&info_t);
    noise = sf_open("noise-sample.wav",SFM_READ,&info_n);
    
    // Errors if opening fails
    if (tango == NULL)
    {
        printf("Failed to open the Tango.\n");
        exit(-1);
    }
    if (noise == NULL)
    {
        printf("Failed to open the noise.\n");
        exit(-1);
    }
    
    // Just equating some things so I don't have to type the structure
    // name every time I want to use these values.
    f_t = info_t.frames;
    sr_t = info_t.samplerate;
    c_t = info_t.channels;
    
    f_n = info_n.frames;
    sr_n = info_n.samplerate;
    c_n = info_n.channels;
    
    num_t_items = f_t*c_t;
    num_n_items = f_n*c_n;
    
    double time = (double)(f_t)/(double)(sr_t);
    
    // Allocate space for the data to be read and read it in.
    buf_t_i = (int *) malloc(num_t_items*sizeof(int));
    num_t = sf_read_int(tango,buf_t_i,num_t_items);
    sf_close(tango);
    
    buf_n_i = (int *) malloc(num_n_items*sizeof(int));
    num_n = sf_read_int(noise,buf_n_i,num_n_items);
    sf_close(noise);
    
    // These are given in int, but I want to use double
    buf_t = (double *) malloc(num_t_items*sizeof(double));
    buf_n = (double *) malloc(num_n_items*sizeof(double));
    
    for(i=0; i<num_t_items; i++) buf_t[i] = (double)(buf_t_i[i]);
    for(i=0; i<num_n_items; i++) buf_n[i] = (double)(buf_n_i[i]);
    
    // We now have the WAV files converted into arrays,
    // which will allow us to find the power spectra.

    // Here i'm taking the Fourier transform of the data and noise.
    gsl_fft_real_radix2_transform(buf_t, 1, num_t_items);
    gsl_fft_real_radix2_transform(buf_n, 1, num_n_items);
    
    // Allocating a bunch of memory.
    double *PSD_t = (double *) malloc(num_t_items/2*sizeof(double));
    double *PSD_n = (double *) malloc(num_n_items/2*sizeof(double));
    double *Gaussian = (double *) malloc(num_n_items/2*sizeof(double));
    double *impulse = (double *) malloc(5000*sizeof(double));
    double *Gaussian_impulse = (double *) malloc(5000*sizeof(double));
    
    // Calculating the spectra of the data and noise.
    for(i=1; i<num_t_items/2; i++) PSD_t[i] = buf_t[i]*buf_t[i] + buf_t[num_t_items-i]*buf_t[num_t_items-i];
    for(i=1; i<num_n_items/2; i++) PSD_n[i] = buf_n[i]*buf_n[i] + buf_n[num_n_items-i]*buf_n[num_n_items-i];
        
    FILE *out;
    FILE *out2;
    
    
    double f;
    
    
    // Opening the files for the spectra and unit impulse smoothing test.
    out = fopen("PSD_T_N.txt","w");
    out2 = fopen("impulse.txt","w");
    
    // Just to make sure I'm not going out of the bounds of my arrays, I'm checking to see if which WAV file is longer.
    int N;
    if(num_t_items<=num_n_items) N = num_t_items;
    if(num_n_items<num_t_items) N = num_n_items;
    
    // Defining my standard deviation for the Gaussian and the normalize Gaussian.
    double sigma = 0.02*sr_t*(double)(N/2);
    for(i=0; i<N/2; i++) Gaussian[i] = (1.0/(sigma*2.507))*exp(-((double)(i)*sr_t)*((double)(i)*sr_t)/(2.0*sigma*sigma));
    
    // Making another standard deviation for the impulse smoothing test.
    sigma = 60.0;
    
    // Calculating the step and Gaussian
    for(i=-2500; i<2500; i++)
    {
        f = (double)(i);
        Gaussian_impulse[i+2500] = (1.0/(sigma*2.507))*exp(-f*f/(2.0*sigma*sigma));
        if(f<-10.0) impulse[i+2500] = 0.0;
        if(f>=-10.0 && f<=10.0) impulse[i+2500] = 1.0/21.0;
        if(f>10.0) impulse[i+2500] = 0.0;
        fprintf(out2, "%.6e %.6e %.6e \n", f, impulse[i+2500], Gaussian_impulse[i+2500]);
    }
    
    // Exporting the spectra file
    for(i=0; i<N/2; i++) fprintf(out, "%.6e %.6e %.6e %.6e \n", (double)(i)*sr_t, PSD_t[i], PSD_n[i], Gaussian[i]);
    fclose(out);
    fclose(out2);
    
    // Allocating more memory
    double *PSD_t_smooth = (double *) malloc(num_t_items/2*sizeof(double));
    double *PSD_n_smooth = (double *) malloc(num_n_items/2*sizeof(double));
    double *impulse_smooth = (double *) malloc(10000*sizeof(double));
    double *Gaussin_long = (double *) malloc(10000*sizeof(double));
    for(i=0; i<10000; i++) Gaussin_long[i] = (1.0/(sigma*2.507))*exp(-(double)(i-2500)*(double)(i-2500)/(2.0*sigma*sigma));
        
    
    // Taking the Fourier transform of the spectra and Gaussian. Since those are all real,
    // I'm using the real transform.
    gsl_fft_real_radix2_transform(PSD_t, 1, N/2);
    gsl_fft_real_radix2_transform(PSD_n, 1, N/2);
    gsl_fft_real_radix2_transform(Gaussian, 1, N/2);
    
    // Multiplying my spectra with my Gaussian.
    for(i=0; i<N/4; i++)
    {
        PSD_t_smooth[i] = PSD_t[i]*Gaussian[i] - PSD_t[N/2-i]*Gaussian[N/2-i];
        PSD_t_smooth[N/2-i] = PSD_t[i]*Gaussian[N/2-i] + PSD_t[N/2-i]*Gaussian[i];
         
        PSD_n_smooth[i] = PSD_n[i]*Gaussian[i] - PSD_n[N/2-i]*Gaussian[N/2-i];
        PSD_n_smooth[N/2-i] = PSD_n[i]*Gaussian[N/2-i] + PSD_n[N/2-i]*Gaussian[i];
    }

    // Using the inverse fft
    gsl_fft_halfcomplex_radix2_inverse(PSD_t_smooth, 1, N/2);
    gsl_fft_halfcomplex_radix2_inverse(PSD_n_smooth, 1, N/2);
    gsl_fft_halfcomplex_radix2_inverse(Gaussian, 1, N/2);
        
    
    // Exporting a file of the smoothed spectra and the Gaussian I used.
    out = fopen("PSD_T_N_smooth_short.txt","w");
    for(i=0; i<N/2; i++) fprintf(out, "%.6e %.6e %.6e %.6e \n", (double)(i)*sr_t, PSD_t_smooth[i], PSD_n_smooth[i], Gaussian[i]);
    fclose(out);
        

    // Calculating the convolution for the impulse with a Gaussian.
    convolve(impulse, 5000, Gaussian_impulse, 5000, impulse_smooth);
    
    // Exporting the convolved impulse.
    out = fopen("impulse_smooth.txt","w");
    for(i=0; i<10000; i++) fprintf(out, "%.6e %.6e %.6e \n", (double)(i)-2500, impulse_smooth[i], Gaussin_long[i]);
    fclose(out);

    
    // Allocating memory, again.
    double *f_filter = (double *) malloc(num_t_items*sizeof(double));
    double *signal = (double *) malloc(num_t_items*sizeof(double));
    double PSD_signal;
    double a;
    
    // Calculating the optimal filter and applying it to my data to get my signal.
    for(i=0; i<N; i++)
    {
        if(i<N/2) a = (PSD_t_smooth[i] - PSD_n_smooth[i])/PSD_t_smooth[i];
        if(i>=N/2) a = (PSD_t_smooth[N-i] - PSD_n_smooth[N-i])/PSD_t_smooth[N-i];
        if(a<0.0) a = 0.0;
        
        f_filter[i] = sqrt(a);
        if(i<N/2) signal[i] = buf_t[i]*f_filter[i];
        if(i>=N/2) signal[i] = buf_t[i]*f_filter[i];
    }
    
    // Exporting a data file with the filter in the frequency domain.
    out = fopen("Filter.txt", "w");
    for(i=0; i<N; i++) fprintf(out, "%.6e %.6e \n", (double)(i)*sr_t, f_filter[i]);
    fclose(out);
    
    // Creating signal spectrum and exporting it.
    out = fopen("PSD_Signal.txt","w");
    for(i=0; i<N/2; i++)
    {
        PSD_signal = signal[i]*signal[i] + signal[N-i]*signal[N-i];
        fprintf(out, "%.6e %.6e \n", (double)(i)*sr_t, PSD_signal);
    }
    fclose(out);
    
    // Taking the inverse fft for the signal.
    gsl_fft_halfcomplex_radix2_inverse(signal, 1, N);
    
    // Exporting the signal.
    out = fopen("signal.txt", "w");
    for(i=1; i<N; i++) fprintf(out, "%.6e %.6e\n", (double)(i)*(time)/((double)(N)), signal[i]);
    fclose(out);
    
    // Taking the real fft of the filter, since it is defined to be real.
    gsl_fft_halfcomplex_radix2_inverse(f_filter, 1, N/2);
    
    // Exporting the filter in the time domain.
    out = fopen("Filter_t.txt", "w");
    for(i=1; i<N/2; i++) fprintf(out, "%.6e %.6e\n", (double)(i)*time/((double)(N/2)), f_filter[i]);
    fclose(out);
    
    
    // Now I want to create my wave file, this code follows from the example code
    // given by the linsndfile library, which is in the folder called make_sine.c
    SNDFILE    *file ;
    SF_INFO    sfinfo ;
    int        k ;
    int    *buffer ;
        
    if (! (buffer = malloc (2 * f_t * sizeof (int))))
    {
        printf ("Malloc failed.\n") ;
        exit (0) ;
    }
        
    memset (&sfinfo, 0, sizeof (sfinfo)) ;
        
    sfinfo.samplerate    = sr_t ;
    sfinfo.frames        = f_t ;
    sfinfo.channels        = 1 ;
    sfinfo.format        = (SF_FORMAT_WAV | SF_FORMAT_PCM_24) ;
        
    if (! (file = sf_open ("signal.wav", SFM_WRITE, &sfinfo)))
    {
        printf ("Error : Not able to open output file.\n") ;
        free (buffer) ;
        return 1 ;
    }
        
    if (sfinfo.channels == 1)
    {
        for (k = 0 ; k < f_t ; k++) buffer [k] = ceil(signal[k]) ;
    }
        
    if (sf_write_int (file, buffer, sfinfo.channels * f_t) != sfinfo.channels * f_t)
    {
        puts (sf_strerror (file)) ;
    }
        
    sf_close(file) ;
    
    // Freeing memory.
    free(buffer) ;
    free(buf_t_i); free(buf_n_i); free(buf_t); free(buf_n); free(PSD_t); free(PSD_n); free(impulse); free(Gaussian_impulse);
    free(PSD_t_smooth); free(PSD_n_smooth); free(impulse_smooth); free(f_filter); free(signal);

}

void convolve(double *Signal, int SignalLen, double *Kernel, int KernelLen, double *Result)
{
    int n;
    int kmin, kmax, k;
    
    for (n = 0; n < SignalLen + KernelLen - 1; n++)
    {
        
        Result[n] = 0;
        
        kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
        kmax = (n < SignalLen - 1) ? n : SignalLen - 1;
        
        for (k = kmin; k <= kmax; k++)
        {
            Result[n] += Signal[k] * Kernel[n - k];
        }
    }
}
