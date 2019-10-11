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
    char errm[50];
    FILE *pf;
    char buff[50], buff2[50], *p;
    int i;

    pf = fopen(parfname, "r");
    if(pf==NULL)
    {
        perror("Error opening file");
        exit(21);
    }
    for(i=0;i<11;i++)
    {
        if(fscanf(pf, "%s", buff)==0)
        {
            perror("Error reading parameter names.");
            exit(3);
        }
        if(strcmp(buff,"datadir") == 0)
        {
            if(fscanf(pf, "%s", pars.datadir)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
        }
        else if(strcmp(buff,"filterfile") == 0)
        {
            if(fscanf(pf, "%s", pars.filterfname)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
        }
        else if(strcmp(buff,"outputdir") == 0)
        {
            if(fscanf(pf, "%s", pars.outputdir)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
        }
        else if(strcmp(buff,"filter_length") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.filterlen = strtol(buff2,&p,10);
        }
        else if(strcmp(buff, "filter_chans") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.filterchans = (int)strtol(buff2, &p, 10);
        }
        else if(strcmp(buff,"amplification") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.ampl = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff,"nchannels") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.nchannels = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff,"firstchan") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.firstchan = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff,"ntiles") == 0)
        {
            if(fscanf(pf, "%s", buff2)==0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.ntiles = (int)strtol(buff2,&p,10);
        }
        else if(strcmp(buff, "tile") == 0)
        {
            if(fscanf(pf, "%s", buff2) == 0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.tile = (int)strtol(buff2, &p, 10);
        }
        else if(strcmp(buff, "pol") == 0)
        {
            if(fscanf(pf, "%s", buff2) == 0)
            {
                sprintf(errm,"Error reading %s\n",buff);
                perror(errm);
                exit(3);
            }
            pars.pol = (int)strtol(buff2, &p, 10);
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
    if(fread(buffer, data_length * sizeof(uint8_t), 1, file)!= 1)
    {
        perror("Error reading data file\n");
        exit(25);
    }

    for(i=0;i<data_length;i++){
        data[i] = buffer[i];
    }
    //free memory
    free(buffer);
}

void actually_read_vcs(FILE *file, uint8_t data[],int data_length, struct parameters pars) /*reads the data from
                                                            the vcs file*/
{
    //open file and variable Initialisation
    uint8_t buffer;
    int i;

    //Read file into buffer as uint8_t
    for(i=0;i<data_length/2;i++){
        if(fread(&buffer, sizeof(uint8_t), 1, file) != 1)
        {
            perror("Error reading data file\n");
            exit(22);
        }
        fseek(file, pars.nchannels*pars.ntiles*2 - 1,SEEK_CUR);
        data[i*2] = buffer & 0xf;
        data[i*2+1] = (buffer >> 4) & 0xf;
    }
}

void read_filter(char *filename, int16_t fdata[],unsigned long filter_length)
//Reads in the data from the specified filter file
{
    //open file and initialise variables
    FILE *fp = fopen(filename,"r");
    float buffer[filter_length];
    char errm[50];
    int i;

    if(fp == NULL)
    {
        sprintf(errm,"The filter file was not found at %s\n",filename);
        perror(errm);
        exit(2);
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
    if(strcmp(pars.outputdir,"notoutfile") == 0)
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
    if(pars.ampl == -1)
    {
        printf("amplification not specified\n");
        a=1;
    }
    if(pars.nchannels == -1)
    {
        printf("nchannels not specified\n");
        a=1;
    }
    if(pars.nchannels == 1)
    {
        printf("The processing is pointless with only one channel...\n");
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
    if(pars.pol == -1)
    {
        printf("polarisation not specified\n");
        a=1;
    }
    if(a)
    {
        perror("Please specify each parameter and try again\n");
        exit(4);
    }
    else
    {
        printf("All parameters specified.\n");
    }
}
