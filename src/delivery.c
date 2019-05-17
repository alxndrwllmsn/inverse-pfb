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
        if(strcmp(buff,"datadir") == 0)
        {
            if(fscanf(pf, "%s", pars.datadir)==0)
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

void read_vcs(FILE *file, uint8_t data[],int data_length) /*reads the data from
                                                            the vcs file*/
{
    //open file and variable Initialisation
    uint8_t *buffer;
    int i;

    //Read file into buffer as uint8_t
    buffer = (uint8_t *)malloc(data_length * sizeof(uint8_t));
    if(fread(buffer, data_length * sizeof(uint8_t), 1, file)==0)
    {
        printf("Error reading data file\n");
        abort();
    }

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

void checkpars(struct parameters pars)
{
    //Initialise variables
    int a = 0;

    //Check each parameter against initialised values
    if(strcmp(pars.datadir,"notdfile") == 0)
    {
        printf("datafile not specified\n");
        a=1;
    }
    if(strcmp(pars.filterfname,"notffile") == 0)
    {
        printf("filterfile not specified\n");
        a=1;
    }
    if(strcmp(pars.outputfname,"notoutfile") == 0)
    {
        printf("outputfile not specified\n");
        a=1;
    }
    if(pars.filterlen == -1)
    {
        printf("filter_length not specified\n");
        a=1;
    }
    if(pars.filterchans == -1)
    {
        printf("filter_chans not specified\n");
        a=1;
    }
    if(pars.nsamples == -1)
    {
        printf("nsamples not specified\n");
        a=1;
    }
    if(pars.nchannels == -1)
    {
        printf("nchannels not specified\n");
        a=1;
    }
    if(pars.firstchan == -1)
    {
        printf("firstchan not specified\n");
        a=1;
    }
    if(pars.ntiles == -1)
    {
        printf("ntiles not specified\n");
        a=1;
    }
    if(pars.tile == -1)
    {
        printf("tile not specified\n");
        a=1;
    }
    if(a)
    {
        printf("Please specify each parameter and try again\n");
        abort();
    }
    else
    {
        printf("All parameters specified.\n");
    }
}

int checkmem(void)
{
    //Initialise variables
    FILE *memfile;
    char buffer[50];
    int memAv, swFree;

    //Read file
    memfile = fopen("/proc/meminfo", "r");
    while(!feof(memfile))
    {
        if(fscanf(memfile, "%s", buffer)==0)
        {
            printf("Error reading memfile.\n");
            abort();
        }
        if(strcmp(buffer,"MemAvailable:") == 0)
        {
            if(fscanf(memfile, "%s", buffer)==0)
            {
                printf("Error reading memory value\n");
                abort();
            }
            memAv = (int)strtol(buffer, NULL, 10);
        }
        else if(strcmp(buffer, "SwapFree:")==0)
        {
            if(fscanf(memfile, "%s", buffer)==0)
            {
                printf("Error reading swap value\n");
                abort();
            }
            swFree = (int)strtol(buffer, NULL, 10);
        }
    }
    fclose(memfile);
    //calculate memory and return
    memAv = memAv + swFree;
    return memAv;
}
