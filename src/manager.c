#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "inttypes.h"
#include "fftw3.h"
#include "delivery.h"
#include "kitchen.h"

int main(int argc, char *argv[])
{
    //Check arguments
    if(argc < 2)
    {
        printf("usage: ipfbrun parameterfilename\n");
        return 1;
    }

    //Initialise variables
    char *parfile = argv[1];
    char fname[100], infotextcat[50], buffer[50], errormessage[100];
    char *infotext = ".info";
    struct parameters pars = {"notdfile","notffile","notoutfile",-1,-1,-1,-1,-1,
                                -1,-1,-1};
    struct dsampled ds;
    int memAv, ntaps, fchans, firstchan, nchans, fact2, i, r, k, n, flip, nsections,
        sectionSize, wholeSection, diff;
    int memUseFactor = 22; // Bytes memory used per bytes input
    long int flength;
    int16_t *fdata;
    uint8_t *chandata;
    int8_t *odata;
    float *data, *rndata, *indata, *predata;
    float tmpr,tmpi,memUse;
    float rmin = 0;
    float imin = 0;
    FILE *test, *test2, *info, *ofile;
    fftw_complex *in, *out, *in1,*out1;
    fftw_plan p, q, m;

    //read parameter file
    pars = getpars(parfile);

    //check parameters
    checkpars(pars);

    printf("datadir: %s\n"
            "filterfile: %s\n"
            "outputdir: %s\n"
            "filterlen: %ld\n"
            "filterchans: %d\n"
            "nsamples: %ld\n"
            "nchannels: %d\n"
            "firstchan: %d\n"
            "ntiles: %d\n"
            "tile: %d\n"
            "pol: %d\n",pars.datadir, pars.filterfname, pars.outputdir,
            pars.filterlen, pars.filterchans, pars.nsamples, pars.nchannels,
            pars.firstchan, pars.ntiles, pars.tile, pars.pol);

    //read filter
    flength = pars.filterlen;
    strcpy(fname, pars.filterfname);

    fdata = (int16_t *)malloc(flength * sizeof *fdata);
    read_filter(fname, fdata, flength);

    // test = fopen("testing/fdatatest.dat", "w");
    // fwrite(fdata, flength * sizeof *fdata, 1, test);
    // fclose(test);

    //import wisdom
    if(exists("ipfbwisdom.ws"))
    {
        if(fftw_import_wisdom_from_filename("ipfbwisdom.ws")==0)
        {
            printf("Wisdom was not loaded correctly, this may take a while\n");
        }
    }
    else
    {
        printf("Wisdom file does not exist, this may take a while. A new\n"
            "wisdom file will be created after this run.\n");
    }

    //determine new sample rate
    fchans = pars.filterchans;
    firstchan = pars.firstchan;
    nchans = pars.nchannels;

    char fnames[nchans][100];
    FILE *dfiles[nchans];

    ds = find_downsampled(fchans,firstchan, nchans);

    //polyphase pad and fft filter
    ntaps = (int)flength/fchans;
    fact2 = ds.factor*2;

    float rdata[fact2],idata[fact2];
    float qrm[fact2][ntaps];
    float rmax[fact2];
    float imax[fact2];

    memset(rmax, 0, fact2 * sizeof(float));
    memset(imax, 0, fact2 * sizeof(float));

    for (r=0;r<fact2;r++)
    {
        for (n=0;n<ntaps;n++)
        {
            qrm[r][n] = (float)fdata[n*fchans + r*((int)fchans/(fact2))];
        }
    }

    // test2 = fopen("testing/qrmtest.dat", "w");
    // printf("factor*2:%d, ntaps:%d\n",fact2,ntaps);
    // fwrite(qrm, fact2*ntaps * sizeof(float), 1, test2);
    // fclose(test2);

    //check flip
    flip = ((int)ds.high/ds.factor+1)%2;

    //check files
    strcpy(infotextcat,pars.datadir);
    strcat(infotextcat,infotext);
    info = fopen(infotextcat,"r");
    if(info == NULL)
    {
        sprintf(errormessage,"Opening %s failed: ",infotextcat);
        perror(errormessage);
        return 1;
    }
    for (i=0;i<nchans;i++)
    {
        if(fscanf(info, "%s", fnames[i])==0)
        {
            perror("Error reading parameter names.\n");
            exit(3);
        }
        sprintf(fname,"%s/%s",pars.datadir,fnames[i]);
        strcpy(fnames[i],fname);
        printf("%s\n",fnames[i]);
    }
    fclose(info);
    //open each data file
    for (i=0;i<nchans;i++)
    {
        dfiles[i] = fopen(fnames[i],"r");
        if(dfiles[i] == NULL)
        {
            sprintf(errormessage,"Could not open file %s\n",fnames[i]);
            perror(errormessage);
            exit(23);
        }
        fseek(dfiles[i], 4096+102400*2*pars.ntiles, SEEK_SET);
        fseek(dfiles[i], 102400*(2*pars.tile+pars.pol), SEEK_CUR);
        printf("Moved marker into position\n");
    }

    //check available memory
    // memUse = fact2*pars.nsamples*memUseFactor;

    //check number of sections (based on memory)
    nsections = 200;
    sectionSize = 51200;
    wholeSection = sectionSize + ntaps;

    //allocate memory for data
    data = calloc(wholeSection * fact2 * 2, sizeof *data);
    chandata = (uint8_t *)malloc(2 * sectionSize *sizeof *chandata);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fact2);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fact2);

    printf("Planning FFT's\n");
    p = fftw_plan_dft_1d(fact2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_free(in);
    fftw_free(out);

    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wholeSection);
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wholeSection);
    diff = 0;

    q = fftw_plan_dft_1d(wholeSection,in1,out1,FFTW_FORWARD, FFTW_ESTIMATE);
    m = fftw_plan_dft_1d(wholeSection,in1,out1,FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_free(in1);
    fftw_free(out1);
    printf("Saving wisdom\n");
    if(fftw_export_wisdom_to_filename("ipfbwisdom.ws")==0)
    {
        printf("Wisdom was not saved correctly.\n");
    }

    rndata = (float *)malloc(wholeSection * sizeof *rndata);
    indata = (float *)malloc(wholeSection * sizeof *indata);
    predata = calloc(2 * ntaps * fact2, sizeof *predata);

    strcpy(infotextcat,pars.outputdir);
    sprintf(buffer,"/out_%d_%d.dat",pars.tile, pars.pol);
    strcat(infotextcat,buffer);
    ofile = fopen(infotextcat,"w");

    odata = (int8_t *)malloc(sectionSize * fact2 * 2 *sizeof *odata);

    //loop over sections-> for each section
    for(i = 0;i<nsections;i++)
    {
        printf("Starting section %d\n",i+1);
        //read section from file
        for (k=0;k<nchans;k++)
        {
            printf("Reading in channel %d\n",k+1);
            // printf("%ld\n",ftell(dfiles[k]));

            read_vcs(dfiles[k], chandata, sectionSize*2);
            for (n=0;n<sectionSize;n++)
            {
                tmpr = (float)(int)chandata[2*n];
                if(tmpr >= 128)
                {
                    tmpr -= 256;
                }
                tmpi = (float)(int)chandata[2*n+1];
                if(tmpi >= 128)
                {
                    tmpi -= 256;
                }
                if(flip)
                {
                    data[(2*(n+ntaps))*fact2+(ds.high - firstchan - k)] = tmpr;
                    data[(2*(n+ntaps))*fact2+(fact2 - (ds.high - firstchan - k))] = tmpr;
                    data[(2*(n+ntaps)+1)*fact2+(ds.high - firstchan - k)] = tmpi;
                    data[(2*(n+ntaps)+1)*fact2+(fact2 - (ds.high - firstchan - k))] = -tmpi;
                }
                else
                {
                    data[(2*(n+ntaps))*fact2+(k+firstchan - ds.low)] = tmpr;
                    data[(2*(n+ntaps))*fact2+(fact2 - (k+firstchan - ds.low))] = tmpr;
                    data[(2*(n+ntaps)+1)*fact2+(k+firstchan - ds.low)] = tmpi;
                    data[(2*(n+ntaps)+1)*fact2+(fact2 - (k+firstchan - ds.low))] = tmpi;
                }

            }
            // if(k==0)
            // {
            //     if (i==0)
            //     {
            //         FILE *test3 = fopen("chandatareadtest.dat", "w");
            //         fwrite(chandata, 2 * sectionSize * sizeof(uint8_t), 1, test3);
            //         fclose(test3);
            //     }
            //     else if(i == 1)
            //     {
            //         FILE *test4 = fopen("chandatareadtest2.dat", "w");
            //         fwrite(chandata, 2 * sectionSize * sizeof(uint8_t), 1, test4);
            //         fclose(test4);
            //     }
            // }
            fseek(dfiles[k], 102400*(2*pars.ntiles-1), SEEK_CUR);

        }

        // FILE *test4 = fopen("testing/datareadtest.dat", "w");
        // fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test4);
        // fclose(test4);
        /*perform ipfb
        {
            perform ifft*/
        printf("Performing iFFT\n");
        for(n=0;n<sectionSize;n++)
        {
            for(k=0;k<fact2;k++)
            {
                rdata[k] = data[(2*(n+ntaps))*fact2+k];
                idata[k] = data[(2*(n+ntaps) + 1)*fact2+k];
            }

            fft(rdata, idata, fact2, rdata, idata, p);
            for(r=0;r<fact2;r++)
            {
                data[((n+ntaps)*2)*fact2 + r] = rdata[r];
                data[((n+ntaps)*2 + 1)*fact2 + r] = idata[r];

            }
        }


        // FILE *test5 = fopen("testing/dataffttest.dat", "w");
        // fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test5);
        // fclose(test5);


            /*prepend extra data unless it is the first section*/
        if (i > 0)
        {
            for(r=0;r<fact2;r++)
            {
                for(n=0;n<ntaps;n++)
                {
                    data[(n*2)*fact2 + r] = predata[(n*2)*fact2 + r];
                    data[(n*2 + 1)*fact2 + r] = predata[(n*2 + 1)*fact2 + r];
                }
            }
        }

        // if (i==1)
        // {
        //     FILE *test6 = fopen("testing/dataprependtest.dat", "w");
        //     fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test6);
        //     fclose(test6);
        // }

            /*perform convolution
            {
                fft section

                multiply with filter

                ifft section
            }*/
        printf("Performing convolution\n");
        for(r=0;r<fact2;r++)
        {
            for(n=0;n<wholeSection;n++)
            {
                rndata[n] = data[(n*2)*fact2 + r];
                indata[n] = data[(n*2 + 1)*fact2 + r];
            }
            for(n=0;n<ntaps;n++)
            {
                predata[(n*2)*fact2 + r] = rndata[n+sectionSize];
                predata[(n*2 + 1)*fact2 + r] = indata[n+sectionSize];
            }
            fftconvolve(rndata, indata, wholeSection, qrm[r], ntaps, rndata, indata, q, m);
            for(n=0;n<sectionSize;n++)
            {
                data[((n+ntaps)*2)*fact2 + r] = rndata[n+ntaps];
                data[((n+ntaps)*2 + 1)*fact2 + r] = indata[n+ntaps];
            }
            rmax[r] = max(rndata, wholeSection);
            rmin = min(rndata, wholeSection);
            if(rmax[r] > -1*(rmin))
            {
                rmax[r] = 128/rmax[r];
            }
            else
            {
                rmax[r] = -127/rmin;
            }
            imax[r] = max(indata, wholeSection);
            imin = min(indata, wholeSection);
            if(imax[r] > -1*(imin))
            {
                imax[r] = 128/imax[r];
            }
            else
            {
                imax[r] = -127/imin;
            }
        }

        // FILE *test7 = fopen("testing/predatatest.dat", "w");
        // fwrite(predata, 2 * ntaps * fact2 * sizeof(float), 1, test7);
        // fclose(test7);
        //
        // FILE *test8 = fopen("testing/dataconvtest.dat", "w");
        // fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test8);
        // fclose(test8);
        //}
        rmin = min(rmax, fact2);
        imin = min(imax, fact2);

        if (rmin < imin)
        {
            imin = rmin;
        }
        for(r=0;r<fact2;r++)
        {
            for(n=0;n<sectionSize;n++)
            {
                odata[(n*fact2 + r)*2] = (int8_t)round(data[((n+ntaps)*2)*fact2 + r]*imin);
                odata[(n*fact2 + r)*2 + 1] = (int8_t)round(data[((n+ntaps)*2 + 1)*fact2 + r]*imin);
            }
        }
        //write section to file
        printf("Writing section %d\n",i+1);
        fwrite(odata, sectionSize * 2 * fact2 *sizeof *odata, 1, ofile);

    }
    //clean up
    printf("Cleaning up\n");
    free(chandata);
    free(rndata);
    free(indata);
    free(fdata);
    free(odata);
    free(predata);
    free(data);

    fftw_destroy_plan(p);
    fftw_destroy_plan(q);
    fftw_destroy_plan(m);
    fftw_cleanup();

    for(int i=0;i<nchans;i++)
    {
        fclose(dfiles[i]);
    }
    fclose(ofile);

}
