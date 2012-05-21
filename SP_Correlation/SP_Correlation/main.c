/**
 * @file main.c
 * @author Bryant Adams
 * @date 3/20/12
 */

#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include "defines.h"
#include <assert.h>
#include <time.h>
#include "tests.h"


/**
 * @brief Entry point for program
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @sideeffect Warns on bad syntax. Passes arguments for processing.
 */
int main (int argc, const char * argv[])
{
    //if(true) runTests();
    //return -1;
    //if(true) makeTestFiles();
    int timestamp = (unsigned)time(NULL);
    int fields;
    
    if(argc%2 != 0 || argc<4)
    {
        printf("Syntax:\n");
        printf("%s {trials} F1a.tdv F2a.tdv F1b.tdv F2b.tdv ...\n", argv[0]);
        return EXIT_FAILURE;
    }
    int trials = atoi(argv[1]);
    if(trials<=0)
    {
        printf("{trials} must be an integer greater than 0\n");
        return EXIT_FAILURE;
    }
    fields = argc-2;

    srand(timestamp); //Seed random number generator

    char fname[100];
    sprintf(fname, "testinfo.%d.report.txt", timestamp);
    
    FILE* output = fopen(fname, "w");

    fprintf(output, "Processing based on command:\n\t");
    for(int i=0; i<argc; i++) fprintf(output, "%s ", argv[i]);
    fprintf(output, "\n\n");

    if(fields%2==0 && fields%3==0) 
        fprintf(output, "WARNING: Unable to infer from number of fields whether you want a Mantel or Partial Mantel test. I'll try both.\n");
    if(fields%2==0)
    {
        fprintf(output, "Processing %d-field Mantel Test on:\n", fields/2);
        for(int i=0; i<fields/2; i++)
        {
            fprintf(output, "\t%s vs %s\n", argv[2*i+2], argv[2*i+3]);
        }
        fprintf(output, "Timestamp: %d\n",timestamp);
        processFilePairs(trials, fields/2, argv, timestamp);
    }
    
    fprintf(output, "\n\n");
    timestamp++;
    /* Partial tests
    if(fields%3==0)
    {
        printf("Processing %d-field Partial Mantel Test on:\n", fields/3);
        for(int i=0; i<fields/3; i++)
        {
            printf("\t%s vs %s given %s\n", argv[3*i+2], argv[3*i+3], argv[3*i+4]);
        }
        printf("Timestamp: %d\n",timestamp);
        processFileTriples(trials, fields/3, argv, timestamp);
    }
    */
    fprintf(output, "\n\n");

    fclose(output);
    printf("\nWriting results to testdata.%d.*\n",timestamp); 

    return EXIT_SUCCESS;
}
    
