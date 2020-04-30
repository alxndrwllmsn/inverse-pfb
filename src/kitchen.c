#include <stdlib.h>
#include "math.h"
#include "inttypes.h"
#include "fftw3.h"
#include "time.h"
#include "kitchen.h"
#include "delivery.h"

float max(float array[], int array_length)
{
    int i;
    float buf = -32000;

    for(i=0;i<array_length;i++){
        if (array[i] > buf){
            buf = array[i];
        }
    }
    return buf;
}

float maxi(int16_t array[], int array_length)
{
    int i;
    float buf = -32000;

    for(i=0;i<array_length;i++){
        if (array[i] > buf){
            buf = array[i];
        }
    }
    return buf;
}

float min(float array[], int array_length)
{
    int i;
    float buf = 32000;

    for(i=0;i < array_length;i++)
    {
        if (array[i] < buf)
        {
            buf = array[i];
        }
    }
    return buf;
}

void convolve(float data[],int nsamples, float filter[], int flength, float outAr[])
{
    //Performs a convolution between a data array and a filter
    int kmin,kmax,k,n;
    int ntotal = nsamples + flength - 1;

    for(n=0;n<ntotal;n++)
    {

        outAr[n] = 0;
        if(n < flength)
        {
            kmin = 0;
        }
        else
        {
            kmin = n - flength + 1;
        }

        if(n >= nsamples)
        {
            kmax = nsamples;
        }
        else
        {
            kmax = n;
        }

        for (k=kmin;k<kmax;k++)
        {
            outAr[n] = outAr[n] + data[k] * filter[n-k];
        }
    }
}

void fftconvolve(float rdata[], float idata[], int nsamples, float filter[], int flength,
    float outAr[], float ioutAr[], fftw_plan p, fftw_plan q)
{
    int ntotal = nsamples;
    int i;

    fftw_complex *ind, *outd, *inf, *outf, *inc, *outc;

    ind = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);
    outd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);

    for(i=0;i<nsamples;i++)
    {
        ind[i][0] = rdata[i];
        ind[i][1] = idata[i];
    }

    fftw_execute_dft(p, ind, outd);
    fftw_free(ind);

    inf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);
    outf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);

    for(i=0;i<flength;i++)
    {
        inf[i][0] = filter[i];
        inf[i][1] = 0;
    }
    for(i=flength;i<ntotal;i++)
    {
        inf[i][0] = 0;
        inf[i][1] = 0;
    }

    fftw_execute_dft(p, inf, outf);
    fftw_free(inf);

    inc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);
    outc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);

    for(i=0;i<ntotal;i++)
    {
        inc[i][0] = outd[i][0]*outf[i][0] - outd[i][1]*outf[i][1];
        inc[i][1] = outd[i][0]*outf[i][1] + outd[i][1]*outf[i][0];
    }

    fftw_free(outd);
    fftw_free(outf);

    fftw_execute_dft(q, inc, outc);
    fftw_free(inc);

    for(i=0;i<ntotal;i++)
    {
        outAr[i] = (float)outc[i][0]/ntotal;
        ioutAr[i] = (float)outc[i][1]/ntotal;
    }

    fftw_free(outc);
}

void rfftconvolve(float rdata[], int nsamples, float filter[], int flength,
    float outAr[], fftw_plan p, fftw_plan q)
{
    int ntotal = nsamples;
    int i;

    fftw_complex *outd, *outf, *inc;
    double *ind, *inf, *outc;

    ind = (double*) fftw_malloc(sizeof(double) * ntotal);
    outd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (int)(ntotal/2+1));

    for(i=0;i<nsamples;i++)
    {
        ind[i] = rdata[i];
    }

    fftw_execute_dft_r2c(p, ind, outd);
    fftw_free(ind);

    inf = (double*) fftw_malloc(sizeof(double) * ntotal);
    outf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (int)(ntotal/2+1));

    for(i=0;i<flength;i++)
    {
        inf[i] = filter[i];
    }
    for(i=flength;i<ntotal;i++)
    {
        inf[i] = 0;
    }

    fftw_execute_dft_r2c(p, inf, outf);
    fftw_free(inf);

    inc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (int)(ntotal/2+1));
    outc = (double*) fftw_malloc(sizeof(double) * ntotal);

    for(i=0;i<(ntotal/2+1);i++)
    {
        inc[i][0] = outd[i][0]*outf[i][0] - outd[i][1]*outf[i][1];
        inc[i][1] = outd[i][0]*outf[i][1] + outd[i][1]*outf[i][0];
    }

    fftw_free(outd);
    fftw_free(outf);

    fftw_execute_dft_c2r(q, inc, outc);
    fftw_free(inc);

    for(i=0;i<ntotal;i++)
    {
        outAr[i] = (float)outc[i]/ntotal;
    }

    fftw_free(outc);
}

void fft(float rdata[], float idata[], int nsamples, float odata[], float oidata[],const fftw_plan p)
{
    fftw_complex *in, *out;
    int i;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nsamples);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nsamples);

    for(i=0;i<nsamples;i++)
    {
        in[i][0] = rdata[i];
        in[i][1] = idata[i];
    }

    fftw_execute_dft(p, in, out);

    for(i=0;i<nsamples;i++)
    {
        odata[i] = out[i][0]/nsamples;
        oidata[i] = out[i][1]/nsamples;
    }

    fftw_free(in);
    fftw_free(out);
}

void fft_c2r(float rdata[], float idata[], int nsamples, float odata[], const fftw_plan p)
{
    fftw_complex *in;
    double *out;
    int i;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nsamples);
    out = (double*) fftw_malloc(sizeof(double) * (2*(nsamples-1)));

    for(i=0;i<nsamples;i++)
    {
        in[i][0] = rdata[i];
        in[i][1] = idata[i];
    }

    fftw_execute_dft_c2r(p, in, out);

    for(i=0;i<(2*nsamples-2);i++)
    {
        odata[i] = out[i]/nsamples;
    }

    fftw_free(in);
    fftw_free(out);
}

struct dsampled find_downsampled(int fchan, int chanloc, int nchans)
{
    struct dsampled ds;

    int ulim = chanloc + nchans - 1;

    int i, half;
    ds.factor = 1;
    ds.low = 0;
    ds.high = (int)fchan/2;
    for(i=0;i<(int)fchan/(2*nchans);i++)
    {
        half = (int)(ds.low+ds.high)/2;
        if(half <= chanloc)
        {
            ds.low = half;
            ds.factor = ds.factor*2;
        }
        else if(half > ulim)
        {
            ds.high = half;
            ds.factor = ds.factor*2;
        }
        else
        {
            break;
        }
    }
    ds.factor = (int)fchan/(2*ds.factor);
    return ds;
}
