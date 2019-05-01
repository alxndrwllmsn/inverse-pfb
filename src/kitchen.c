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
    int ntotal = nsamples + flength - 1;
    int hflength = (int)flength/2;
    int hntotal = (int)ntotal/2;
    int i;

    fftw_complex *ind, *outd, *inf, *outf, *inc, *outc;

    ind = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);
    outd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntotal);

    for(i=0;i<nsamples;i++)
    {
        ind[i][0] = rdata[i];
        ind[i][1] = idata[i];
    }
    for(i=nsamples;i<ntotal;i++)
    {
        ind[i][0] = 0;
        ind[i][1] = 0;
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

    for(i=hflength;i<hntotal;i++)
    {
        outAr[i-hflength] = (float)outc[i + hntotal][0]/ntotal;
        ioutAr[i-hflength] = (float)outc[i + hntotal][1]/ntotal;
    }
    for(i=hntotal;i<hflength + nsamples;i++)
    {
        outAr[i - hflength] = (float)outc[i - hntotal][0]/ntotal;
        ioutAr[i-hflength] = (float)outc[i - hntotal][1]/ntotal;
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
        odata[i] = out[i][0]/sqrt(nsamples);
        oidata[i] = out[i][1]/sqrt(nsamples);
    }

    fftw_free(in);
    fftw_free(out);
}

void inverse_pfb(uint8_t data[], int16_t filter[], int nsamples, int filter_length,
    int fchannels, int nchannels, int8_t outData[], int chanloc, struct dsampled ds)
{
    //initialise variables
    int n,k,r, flip;
    float *wdata;     //[ds.factor*2][2][nsamples]
    int ntaps = (int)filter_length/fchannels;
    int fact2 = ds.factor*2;
    float qrm[fact2][ntaps];
    float rdata[fact2];
    float idata[fact2];
    float *rndata, *indata;
    float tmpr, tmpi;
    float rmax[fact2], imax[fact2];
    float rmin = 0;
    float imin = 0;

    clock_t start, diff;
    int msec;

    //allocate array
    wdata = calloc(fact2*2*nsamples, sizeof *wdata);
    rndata = (float *)malloc(nsamples * sizeof *rndata);
    indata = (float *)malloc(nsamples * sizeof *indata);

    //polyphase structure
        //filter
    for(r=0;r<fact2;r++)
    {
        for(n=0;n<ntaps;n++)
        {
            qrm[r][n] = (float) filter[n*fchannels + r*((int)fchannels/(fact2))];
        }
    }
        //determine whether to flip or not
    flip = ((int)ds.high/ds.factor)%2+1;
        //data
    for(n=0;n<nsamples;n++)
    {
        for(k=0;k<nchannels;k++)
        {
            tmpr = (float)(int)data[(nchannels*n+k)*2];
            if(tmpr >= 128)
            {
                tmpr = tmpr - 256;
            }
            tmpi = (float)(int)data[(nchannels*n+k)*2+1];
            if(tmpi >= 128)
            {
                tmpi = tmpi - 256;
            }
            if(flip)
            {
                wdata[(2*n)*fact2+(ds.high - chanloc - k)] = tmpr;
                wdata[(2*n)*fact2+(fact2 - (ds.high - chanloc - k))] = tmpr;
                wdata[(2*n+1)*fact2+(ds.high - chanloc - k)] = tmpi;
                wdata[(2*n+1)*fact2+(fact2 - (ds.high - chanloc - k))] = -tmpi;
            }
            else
            {
                wdata[(2*n)*fact2+(k+chanloc - ds.low)] = tmpr;
                wdata[(2*n)*fact2+(fact2 - (k+chanloc - ds.low))] = tmpr;
                wdata[(2*n+1)*fact2+(k+chanloc - ds.low)] = tmpi;
                wdata[(2*n+1)*fact2+(fact2 - (k+chanloc - ds.low))] = tmpi;
            }
        }
    }

    start = clock();

    //ifft in terms of channel and branch
    fftw_complex *in,*out;
    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fact2);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fact2);

    //check for wisdom
    int ws = exists("wisdom.ws");
    if(ws)
    {
        fftw_import_wisdom_from_filename("wisdom.ws");
    }

    for(n=0;n<nsamples;n++)
    {
        for(k=0;k<fact2;k++)
        {
            rdata[k] = wdata[(2*n)*fact2+k];
            idata[k] = wdata[(2*n + 1)*fact2+k];
        }
        if(n==0)
        {
            p = fftw_plan_dft_1d(fact2, in, out, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
        }
        fft(rdata, idata, fact2, rdata, idata, p);
        for(r=0;r<fact2;r++)
        {
            wdata[(n*2)*fact2 + r] = rdata[r];
            wdata[(n*2 + 1)*fact2 + r] = idata[r];
        }
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("fft data: %d s %d ms\n",msec/1000,msec%1000);
    start = clock();

    //convolve with filter
    fftw_complex *in1,*out1;
    fftw_plan q, m;
    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nsamples);
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nsamples);
    diff = 0;

    q = fftw_plan_dft_1d(nsamples,in1,out1,FFTW_FORWARD, FFTW_EXHAUSTIVE);
    m = fftw_plan_dft_1d(nsamples,in1,out1,FFTW_BACKWARD, FFTW_EXHAUSTIVE);

    //export wisdom
    if(ws == 0)
    {
        fftw_export_wisdom_to_filename("wisdom.ws");
    }

    for(r=0;r<fact2;r++)
    {
        for(n=0;n<nsamples;n++)
        {
            rndata[n] = wdata[(n*2)*fact2 + r];
            indata[n] = wdata[(n*2 + 1)*fact2 + r];
        }
        fftconvolve(rndata, indata, nsamples, qrm[r], ntaps, rndata, indata, q, m);
        rmax[r] = max(rndata, nsamples);
        rmin = min(rndata, nsamples);
        if(rmax[r] > -1*(rmin))
        {
            rmax[r] = 128/rmax[r];
        }
        else
        {
            rmax[r] = -127/rmin;
        }
        imax[r] = max(indata, nsamples);
        imin = min(indata, nsamples);
        if(imax[r] > -1*(imin))
        {
            imax[r] = 128/imax[r];
        }
        else
        {
            imax[r] = -127/imin;
        }

        for(n=0;n<nsamples;n++)
        {
            wdata[(n*2)*fact2 + r] = rndata[n];
            wdata[(n*2 + 1)*fact2 + r] = indata[n];
        }
    }


    rmin = min(rmax, fact2);
    imin = min(imax, fact2);

    if (rmin < imin)
    {
        imin = rmin;
    }
    for(r=0;r<fact2;r++)
    {
        for(n=0;n<nsamples;n++)
        {
            outData[(n*fact2 + r)*2] = (int8_t)round(wdata[(n*2)*fact2 + r]*imin);
            outData[(n*fact2 + r)*2 + 1] = (int8_t)round(wdata[(n*2 + 1)*fact2 + r]*imin);
        }
    }

    free(rndata);
    free(indata);

    fftw_destroy_plan(q);
    fftw_destroy_plan(m);
    fftw_free(in1);
    fftw_free(out1);
    fftw_cleanup();

    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("convolve data and filter: %d s %d ms\n",msec/1000,msec%1000);
    start = clock();

    free(wdata);

}

struct dsampled find_downsampled(int fchan, int chanloc, int nchans)
{
    struct dsampled ds;

    int ulim = chanloc + nchans;

    int i, half;
    ds.factor = 1;
    ds.low = 0;
    ds.high = (int)fchan/2;
    for(i=0;i<(int)fchan/(2*nchans);i++)
    {
        half = (int)(ds.low+ds.high)/2;
        if(half < chanloc)
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
