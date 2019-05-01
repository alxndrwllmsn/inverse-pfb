#ifndef _DELIVERY_H_
#define _DELIVERY_H_

struct parameters
{
    char datafname[50];
    char filterfname[50];
    char outputfname[50];
    long filterlen;
    int filterchans;
    long nsamples;
    int nchannels;
    int firstchan;
    int ntiles;
    int tile;
};

struct parameters getpars(char *parfname);
int read_vcs(char *filename, uint8_t data[], int data_length);
int read_filter(char *filename, int16_t fdata[], unsigned long filter_length);
int exists(const char *fname);
int write_output(char filename[], int8_t array[], int arraysize, char mode[2]);

#endif
