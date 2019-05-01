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
    struct parameters pars;

    //read parameter file
    // pars = getpars(parfile);

    //check parameters
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



    //check number of sections (based on memory)

    //read filter

    //import wisdom

    //polyphase pad and fft filter

    //determine new sample rate

    //open file

    //loop over sections-> for each section
        //read section from file

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

    //clean up
}
