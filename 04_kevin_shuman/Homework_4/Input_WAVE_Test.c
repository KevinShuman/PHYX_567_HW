// This code comes from the following website:
// https://ubuntuforums.org/showthread.php?t=968690

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <string.h>

// Compile with
// gcc -o Input_WAVE_Test Input_WAVE_Test.c -lgsl -lgslcblas -lm -lsndfile

int main()
{
    SNDFILE *sf;
    SF_INFO info;
    int num_channels;
    int num, num_items;
    int *buf;
    int f,sr,c;
    int i,j;
    FILE *out;
    
    /* Open the WAV file. */
    info.format = 0;
    sf = sf_open("Lacumparsita-noisy.wav",SFM_READ,&info);
    if (sf == NULL)
    {
        printf("Failed to open the file.\n");
        exit(-1);
    }
    
    /* Print some of the info, and figure out how much data to read. */
    f = info.frames;
    sr = info.samplerate;
    c = info.channels;
    
    printf("frames=%d\n",f);
    printf("samplerate=%d\n",sr);
    printf("channels=%d\n",c);
    
    num_items = f*c;
    printf("num_items=%d\n",num_items);
    
    
    /* Allocate space for the data to be read, then read it. */
    buf = (int *) malloc(num_items*sizeof(int));
    num = sf_read_int(sf,buf,num_items);
    sf_close(sf);
    
    
    printf("Read %d items\n",num);
    
    
    /* Write the data to filedata.txt. */
    out = fopen("filedata.txt","w");
    for (i = 0; i < num; i += c)
    {
        fprintf(out, "%f ", (double)(i)/(double)(sr));
        for (j = 0; j < c; ++j) fprintf(out,"%d ",buf[i+j]);
        fprintf(out,"\n");
    }
    
    SNDFILE    *file ;
    SF_INFO    sfinfo ;
    int        k ;
    int    *buffer ;
    
    if (! (buffer = malloc (2 * f * sizeof (int))))
    {
        printf ("Malloc failed.\n") ;
        exit (0) ;
    } ;
    
    memset (&sfinfo, 0, sizeof (sfinfo)) ;
    
    sfinfo.samplerate    = sr ;
    sfinfo.frames        = f ;
    sfinfo.channels        = 1 ;
    sfinfo.format        = (SF_FORMAT_WAV | SF_FORMAT_PCM_24) ;
    
    if (! (file = sf_open ("signal_1.wav", SFM_WRITE, &sfinfo)))
    {    printf ("Error : Not able to open output file.\n") ;
        free (buffer) ;
        return 1 ;
    } ;
    
    if (sfinfo.channels == 1)
    {    for (k = 0 ; k < f ; k++)
        buffer [k] = buf[k] ;
    }
    
    if (sf_write_int (file, buffer, sfinfo.channels * f) !=
        sfinfo.channels * f)
        puts (sf_strerror (file)) ;
    
    sf_close (file) ;
    free (buffer) ;
    
    fclose(out);
    return 0;
    }
