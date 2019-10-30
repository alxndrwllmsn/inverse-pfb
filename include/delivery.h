#ifndef _DELIVERY_H_
#define _DELIVERY_H_

struct parameters
{
    char datadir[100];
    char filterfname[100];
    char outputdir[100];
    long filterlen;
    int filterchans;
    int ampl;
    int nchannels;
    int firstchan;
    int ntiles;
    int tile;
    int pol;
};

struct parameters getpars(char *parfname);
void read_vcs(FILE *file, uint8_t data[], int data_length);
void actually_read_vcs(FILE *file, uint8_t data[],int data_length, struct parameters pars);
void read_filter(char *filename, int16_t fdata[], unsigned long filter_length);
int exists(const char *fname);
void checkpars(struct parameters pars);
void write_output(char filename[], int8_t array[], int arraysize, char mode[2]);
#endif
