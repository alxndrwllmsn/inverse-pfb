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
    if(argc < 3)
    {
        printf("usage: ipfbrun parameterfilename vcs\n");
        return 1;
    }

    //Initialise variables
    char *parfile = argv[1];
    int vcs = strtol(argv[2],NULL,10);
    char fname[200], infotextcat[100], buffer[100], errormessage[100];
    char *infotext = ".info";
    struct parameters pars = {"notdfile","notffile","notoutfile",-1,-1,-1,-1,-1,
                                -1,-1,-1};
    struct dsampled ds;
    int ntaps, fchans, firstchan, nchans, fact2, i, r, k, n, flip, nsections,
        sectionSize, wholeSection;
    long int flength;
    int16_t *fdata;
    uint8_t *chandata;
    int8_t *odata;
    float *data, *rndata, *indata, *predata;
    float tmpr,tmpi,fmaxi;
    int imin = 0;
    FILE *info, *ofile;
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
            "amplification: %d\n"
            "nchannels: %d\n"
            "firstchan: %d\n"
            "ntiles: %d\n"
            "tile: %d\n"
            "pol: %d\n",pars.datadir, pars.filterfname, pars.outputdir,
            pars.filterlen, pars.filterchans, pars.ampl, pars.nchannels,
            pars.firstchan, pars.ntiles, pars.tile, pars.pol);

    //read filter
    flength = pars.filterlen;
    strcpy(fname, pars.filterfname);

    fdata = (int16_t *)malloc(flength * sizeof *fdata);
    read_filter(fname, fdata, flength);
    fmaxi = maxi(fdata,flength);
    #ifdef DEBUG
    {
        FILE *test = fopen("fdatatest.dat", "w");
        fwrite(fdata, flength * sizeof *fdata, 1, test);
        fclose(test);
    }
    #endif

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
    printf("Downsampling to %d channels (with conjugate),\nLow channel: %d, High channel: %d\n",ds.factor*2, ds.low, ds.high );
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

    #ifdef DEBUG
    {
        FILE *test2 = fopen("qrmtest.dat", "w");
        printf("factor*2:%d, ntaps:%d\n",fact2,ntaps);
        fwrite(qrm, fact2*ntaps * sizeof(float), 1, test2);
        fclose(test2);
    }
    #endif

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
        if(vcs == 1)
        {
            fseek(dfiles[i], ((i*pars.ntiles + pars.tile)*2 + pars.pol)*2, SEEK_SET);
            printf("%d\n",ftell(dfiles[i]));
        }
        else
        {
            fseek(dfiles[i], 4096+102400*2*pars.ntiles, SEEK_SET);
            fseek(dfiles[i], 102400*(2*pars.tile+pars.pol), SEEK_CUR);
        }
        printf("Moved marker into position\n");
    }

    //check number of sections (based on memory)
    if(vcs== 1)
    {
        nsections = 25;
        sectionSize = 400;
    }
    else
    {
        nsections = 25; //Change this back to 200 after testing
        sectionSize = 51200;
    }
    wholeSection = sectionSize + ntaps;

    //allocate memory for data
    data = calloc(wholeSection * fact2 * 2, sizeof *data);
    chandata = (uint8_t *)malloc(2 * sectionSize *sizeof *chandata);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fact2);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fact2);

    printf("Planning FFT's\n");
    p = fftw_plan_dft_1d(fact2, in, out, FFTW_BACKWARD, FFTW_EXHAUSTIVE);//*This needs to change

    fftw_free(in);
    fftw_free(out);

    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wholeSection);
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wholeSection);

    q = fftw_plan_dft_1d(wholeSection,in1,out1,FFTW_FORWARD, FFTW_EXHAUSTIVE);
    m = fftw_plan_dft_1d(wholeSection,in1,out1,FFTW_BACKWARD, FFTW_EXHAUSTIVE);

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
    printf("%s\n",infotextcat);
    ofile = fopen(infotextcat,"w");
    // sprintf(buffer,"%s/norms_%d_%d.txt",pars.outputdir,pars.tile,pars.pol);
    // norms = fopen(buffer,"w");

    odata = (int8_t *)malloc(sectionSize * fact2 * 2 *sizeof *odata);


    //loop over sections-> for each section
    for(i = 0;i<nsections;i++)
    {
        printf("Starting section %d\n",i+1);
        //read section from file
        for (k=0;k<nchans;k++)
        {
            printf("Reading in channel %d\n",k+1);
            if(vcs==1)
            {
                actually_read_vcs(dfiles[k], chandata, sectionSize*2, pars);
            }
            else
            {
                read_vcs(dfiles[k], chandata, sectionSize*2);
            }
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
                    data[(2*(n+ntaps))*fact2+(ds.high - firstchan - k-1)] = tmpr;
                    data[(2*(n+ntaps)+1)*fact2+(ds.high - firstchan - k-1)] = tmpi;

                    if (k < (ds.high-ds.low-1))
                    {
                        data[(2*(n+ntaps))*fact2+(fact2 - (ds.high - firstchan - k-1))] = tmpr;
                        data[(2*(n+ntaps)+1)*fact2+(fact2 - (ds.high - firstchan - k-1))] = -tmpi;
                    }
                }
                else
                {
                    data[(2*(n+ntaps))*fact2+(k+firstchan - ds.low)] = tmpr;
                    data[(2*(n+ntaps)+1)*fact2+(k+firstchan - ds.low)] = tmpi;
                    if (k > 0)
                    {
                        data[(2*(n+ntaps))*fact2+(fact2 - (k+firstchan - ds.low))] = tmpr;
                        data[(2*(n+ntaps)+1)*fact2+(fact2 - (k+firstchan - ds.low))] = -tmpi;
                        int temp = (2*(n+ntaps)+1)*fact2+(fact2 - (k+firstchan - ds.low));
                        if (temp >= wholeSection * fact2 * 2)
                        {
                            printf("index:%d n:%d k:%d i:%d\n",temp,n,k,i);
                        }
                    }
                }


            }
            #ifdef DEBUG
            {
                if((k==0))
                {
                    FILE *test3 = fopen("chandatareadtest.dat", "w");
                    fwrite(chandata, 2 * sectionSize * sizeof(uint8_t), 1, test3);
                    fclose(test3);
                }
            }
            #endif
            if (vcs!=1)
            {
                fseek(dfiles[k], 102400*(2*pars.ntiles-1), SEEK_CUR);
            }
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
        #ifdef DEBUG
        {
            if (i == 0)
            {
                FILE *test5 = fopen("dataffttest.dat", "w");
                fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test5);
                fclose(test5);
            }
        }
        #endif


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
        #ifdef DEBUG
        {
            if (i==1)
            {
                FILE *test6 = fopen("dataprependtest.dat", "w");
                fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test6);
                fclose(test6);
            }
        }
        #endif

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
                data[((n+ntaps)*2)*fact2 + r] = rndata[n+ntaps]*pars.ampl/fmaxi;
                data[((n+ntaps)*2 + 1)*fact2 + r] = indata[n+ntaps]*pars.ampl/fmaxi;
            }
        }

        #ifdef DEBUG
        {
            FILE *test7 = fopen("predatatest.dat", "w");
            fwrite(predata, 2 * ntaps * fact2 * sizeof(float), 1, test7);
            fclose(test7);

            FILE *test8 = fopen("dataconvtest.dat", "w");
            fwrite(data, 2 * wholeSection * fact2 * sizeof(float), 1, test8);
            fclose(test8);
        }
        #endif
        for(r=0;r<fact2;r++)
        {
            for(n=0;n<sectionSize;n++)
            {
                if(fabs(data[((n+ntaps)*2)*fact2 + r]) > 127)
                {
                    imin++;
                    printf("Over %f\n",fabs(data[((n+ntaps)*2)*fact2 + r]));
                    exit(100);
                }
                odata[(n*fact2 + r)*2] = (int8_t)round(data[((n+ntaps)*2)*fact2 + r]);
                odata[(n*fact2 + r)*2 + 1] = (int8_t)round(data[((n+ntaps)*2 + 1)*fact2 + r]);
            }
        }
        memset(data,0,2 * wholeSection * fact2 * sizeof *data);
        //write normalisation factor to file
        // fprintf(norms, "%d\n",imin);
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
    // fclose(norms);

}
