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
        return -1;
    }

    //Initialise variables
    char *parfile = argv[1];
    struct parameters pars = {"notdfile","notffile","notoutfile",-1,-1,-1,-1,-1,-1,-1};
    struct dsampled ds;
    int memAv, ntaps;
    int16_t *fdata;
    uint8_t *data;

    //read parameter file
    pars = getpars(parfile);

    //check parameters
    checkpars(pars);

    printf("datafile: %s\n"
            "filterfile: %s\n"
            "outputfile: %s\n"
            "filterlen: %ld\n"
            "filterchans: %d\n"
            "nsamples: %ld\n"
            "nchannels: %d\n"
            "firstchan: %d\n"
            "ntiles: %d\n"
            "tile: %d\n",pars.datafname, pars.filterfname, pars.outputfname,
            pars.filterlen, pars.filterchans, pars.nsamples, pars.nchannels,
            pars.firstchan, pars.ntiles, pars.tile);

    //read filter
    long int flength = pars.filterlen;
    char *fname = pars.filterfname;

    fdata = (int16_t *)malloc(flength * sizeof *fdata);
    read_filter(fname, fdata, flength);

    FILE *test = fopen("testing/ftest.dat", "w");
    fwrite(fdata, flength * sizeof *fdata, 1, test);
    fclose(test);

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
    int fchans = pars.filterchans;
    int firstchan = pars.firstchan;
    int nchans = pars.nchannels;

    ds = find_downsampled(fchans,firstchan, nchans);

    //polyphase pad and fft filter
    ntaps = (int)flength/fchans;
    int fact2 = ds.factor*2;

    float qrm[fact2][ntaps];

    for (int r=0;r<fact2;r++)
    {
        for (int n=0;n<ntaps;n++)
        {
            qrm[r][n] = (float)fdata[n*fchans + r*((int)fchans/(fact2))];
        }
    }

    FILE *test2 = fopen("testing/ftest2.dat", "w");
    printf("factor*2:%d, ntaps:%d\n",fact2,ntaps);
    fwrite(qrm, fact2*ntaps * sizeof(float), 1, test2);
    fclose(test2);

    //check flip
    int flip = ((int)ds.high/ds.factor+1)%2;

    //open file
    FILE *dfile = fopen(pars.datafname,"r");
    if (dfile==NULL)
    {
        printf("The datafile was not found\n");
        abort();
    }

    //check available memory
    memAv = checkmem();
    printf("Memory and Swap available: %d kB\n", memAv);

    //check number of sections (based on memory)
    int nsections = 1;
    int sectionSize = (int)pars.nsamples/nsections;

    //allocate memory for data
    data = (uint8_t *)malloc(sectionSize * nchannels *sizeof *data);
    //loop over sections-> for each section
    for(int i = 0;i<nsections;i++)
    {
        //read section from file
        read_vcs(dfile, data, sectionSize);

        //flip or not

        /*perform ipfb
        {
            perform ifft

            prepend extra data (zeros if first section)

            perform convolution
            {
                fft section

                multiply with filter

                ifft section
            }

            remove prepended info

            save end samples for next section
        }*/

        //write section to file

    }
    //clean up
    fclose(dfile);
    free(fdata);
}
