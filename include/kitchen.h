#ifndef _KITCHEN_H_
#define _KITCHEN_H_

struct dsampled
{
    int factor;
    int low;
    int high;
};

float max(float array[], int array_length);
float min(float array[], int array_length);
void convolve(float data[],int nsamples, float filter[], int flength,
    float outAr[]);
void fftconvolve(float rdata[], float idata[], int nsamples, float filter[],
    int flength,float outAr[], float ioutAr[], fftw_plan p, fftw_plan q);
void fft(float rdata[], float idata[], int nsamples, float odata[],
    float oidata[], const fftw_plan p);
struct dsampled find_downsampled(int fchan, int chanloc, int nchans);

#endif
