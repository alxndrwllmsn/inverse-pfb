#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "inttypes.h"
#include "delivery.h"

struct parameters getpars(char *parfname) /*returns the values from a parameter
                                            file*/
{
    //Variable Initialisation
    struct parameters pars;
    FILE *pf;
    char buff[50], buff2[50], *p;
    int i;

    pf = fopen(parfname, "r");
    for(i=0;i<10;i++)
    {
        if(fscanf(pf, "%s", buff)==0)
        {
            printf("Error reading parameter names.\n");
            abort();
        }
        if(strcmp(buff,"datafile") == 0)
        {
            if(fscanf(pf, "%s", pars.datafname)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
        }
        else if(strcmp(buff,"filterfile") == 0)
        {
            if(fscanf(pf, "%s", pars.filterfname)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
        }
        else if(strcmp(buff,"outputfile") == 0)
        {
            if(fscanf(pf, "%s", pars.outputfname)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
        }
        else if(strcmp(buff,"filter_length") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.filterlen = strtol(buff2,&p,10);
        }
        else if(strcmp(buff, "filter_chans") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.filterchans = (int)strtol(buff2, &p, 10);
        }
        else if(strcmp(buff,"nsamples") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.nsamples = strtol(buff2,&p,10);
        }
        else if(strcmp(buff,"nchannels") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.nchannels = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff,"firstchan") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.firstchan = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff,"ntiles") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.ntiles = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff, "tile") == 0)
        {
            if(fscanf(pf, "%s", buff2) == 0)
            {
                printf("Error reading %s\n",buff );
                abort();
            }
            pars.tile = (int)strtol(buff2, &p, 10);
        }
    }
    fclose(pf);
    return pars;
}

void read_vcs(char *filename, uint8_t data[],int data_length) /*reads the data from
                                                            the vcs file*/
{
    //open file and variable Initialisation
    FILE *fp = fopen(filename,"r");
    uint8_t *buffer;
    unsigned long fileLen;
    int i;

    //check length of file and compare to size of array
    fseek(fp, 0, SEEK_END);
    fileLen = ftell(fp);
    if(fileLen != data_length){
        printf("The entered dimensions do not match the number of datapoints "
        "within the file. Please check that the dimensions equal the file "
        "length %ld. The dimensions entered have %d datapoints.\n",fileLen,
        data_length);
        abort();
    }
    fseek(fp, 0, SEEK_SET);

    //Read file into buffer as uint8_t
    buffer = (uint8_t *)malloc(fileLen * sizeof(uint8_t));
    if(fread(buffer, fileLen * sizeof(uint8_t), 1, fp)==0)
    {
        printf("Error reading data file\n");
        abort();
    }
    fclose(fp);

    //convert uint8_t to int and save to data (not the most efficients)
    for(i=0;i<data_length;i++){
        data[i] = buffer[i];
    }
    //free memory
    free(buffer);
}

void read_filter(char *filename, int16_t fdata[],unsigned long filter_length)
//Reads in the data from the specified filter file
{
    //open file and initialise variables
    FILE *fp = fopen(filename,"r");
    float buffer[filter_length];
    int i;

    if(fp == NULL)
    {
        printf("The filter file was not found at %s\n",filename );
        abort();
    }

    //read file and save to fdata
    for(i=0;i<filter_length;i++){
        if(fscanf(fp,"%f", &buffer[i]) == EOF) break;
        fdata[i]=(int16_t) buffer[i];
    }

    //close file
    fclose(fp);
}

int exists(const char *fname)
{
    FILE *file;
    if ((file = fopen(fname, "r")))
    {
        fclose(file);
        return 1;
    }
    return 0;
}

void write_output(char filename[],int8_t array[],int arraysize,char mode[2])
{
    //initialise variables
    FILE *of;

    //open file
    of = fopen(filename, mode);

    //write
    printf("number of elements written: %d\n",arraysize );
    fwrite(array,arraysize * sizeof *array, 1, of);
    fclose(of);
}
