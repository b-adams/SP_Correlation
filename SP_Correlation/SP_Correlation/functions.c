/**
 * @file functions.c
 * @author Bryant Adams
 * @date 3/20/12
 */


#include <stdio.h>
#include "functions.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define VERBOSE false

#pragma mark Utility

float* allocateArrayOfFloats(int howmany)
{
    return malloc(howmany * sizeof(float));
}

float** allocateArrayOfArraysOfFloats(int howmany)
{
    return malloc(howmany * sizeof(float*));
}

int* allocateArrayOfInts(int howmany)
{
    return malloc(howmany * sizeof(int));
}

Field** allocateArrayOfFields(int howmany);
Field** allocateArrayOfFields(int howmany)
{
    return malloc(howmany * sizeof(Field));
}

int randInRange(int lo, int hi)
{
    assert(lo<=hi);
    
    int range = hi-lo+1;
    int num = rand()%range;
    num = num + lo;
    return num;
}
bool inOrder(float l, float r)
{
    return l >= r;
}
void swapF(float* l, float* r)
{
    assert(l!=NULL);
    assert(r!=NULL);

    float temp;
    temp = *l;
    *l = *r;
    *r = temp;
}
void swapI(int* l, int* r)
{
    assert(l!=NULL);
    assert(r!=NULL);

    int temp;
    temp = *l;
    *l = *r;
    *r = temp;
}


#pragma mark Permutation

////////////////////////////////////////////////////
// Permtation adapted from Fisher-Yates Shuffle
// http://en.wikipedia.org/wiki/Fisher–Yates_shuffle
Perm* allocatePermutation(void)
{
    return malloc(sizeof(Perm));
}

Perm* makePerm(int size, int seed)
{
    assert(size>0);

    if(seed>0) srand(seed);
    
    Perm* thePerm = allocatePermutation();
    thePerm->size = size;
    thePerm->index = allocateArrayOfInts(size);
    
    int j;

    thePerm->index[0] = 0;
    for(int i=1; i<size; i++)
    {
        if(!seed)
        {
            thePerm->index[i] = i;
        } else {
            j = randInRange(0, i);
            thePerm->index[i] = thePerm->index[j];
            thePerm->index[j] = i;
        } 
    }
    return thePerm;
}

void modifyPermPermutify(Perm* thePerm, int seed)
{
    assert(thePerm!=NULL);

    if(seed>0) srand(seed);
    int j;
        
    //Including 0 because the last swap hurts nobody,
    //and with no seed the last step finishes setting the
    //permutation to the identity
    for(int i=thePerm->size-1; i>=0; i--)
    {
        if(!seed)
        {
            thePerm->index[i] = i;
        } else {
            j = randInRange(0, i);
            swapI(&thePerm->index[i], &thePerm->index[j]);
        } 
    }
}

#pragma mark Fields

Field* allocateField(void)
{
    return malloc(sizeof(Field));
}

Field* makeRandomField(int samples)
{
    assert(samples>0);

    Field* theField = allocateField();
    theField->samples = samples;
    theField->fieldnum = randInRange(90, 99);
    theField->element = allocateArrayOfArraysOfFloats(samples);
    theField->hasFlatVersion = false;
    theField->perm = makePerm(samples, SEED_IDENTITY);
//    theField->isCentered = false;
//    theField->isRanked = false;
//    theField->isRankBased = false;
//    theField->isRaw = true;
    
    for(int x=0; x<samples; x++)
    {
        theField->element[x] = allocateArrayOfFloats(samples);
    }
    for(int x=0; x<samples; x++)
    {
        for(int y=x; y<samples; y++)
        {
            if(x==y) theField->element[x][y] = 0;
            else
            {
                //theField->element[x][y] = randInRange(0, 99)/10.0;
                theField->element[x][y] = randInRange(0, 9)/1.0;
                theField->element[y][x] = theField->element[x][y];
            }
        }
    }
    return theField;
}

Field* makeFieldFromTDV(const char* filename)
{
    assert(filename!=NULL);

    int fieldnum, samples, scan;
    FILE *theFile = fopen(filename, "r");
    assert(theFile); //Gotta have a file to process
    scan=fscanf(theFile, "<Field %d:Samples %d>",
                &fieldnum, &samples);
    if(scan!=2)
    {
        printf("ERROR: Bad input file <%s> missing header.\n", filename);
        fclose(theFile);
    } else {
        printf("Reading [%s] field=%d, samples=%d\n", filename, fieldnum, samples);
    }
    assert(scan==2);
    
    Field* theField = allocateField();
    theField->samples = samples;
    theField->fieldnum = fieldnum;
    theField->element = allocateArrayOfArraysOfFloats(samples);
    theField->hasFlatVersion = false;
    theField->perm = makePerm(samples, SEED_IDENTITY);
//    theField->isCentered = false;
//    theField->isRanked = false;
//    theField->isRankBased = false;
//    theField->isRaw = true;
    
    for(int x=0; x<samples; x++)
    {
        theField->element[x] = allocateArrayOfFloats(samples);
        for(int y=0; y<samples; y++)
        {
            fscanf(theFile, "%f", &theField->element[x][y]);
            //Funny message when main diagonal nonzero
            if(x==y && theField->element[x][y]!=0) printf("HEY! Nonzero on main diagonal.\n");
        }
    }
    fclose(theFile);
    return theField;
}

void saveFieldToTDV(const char* filename, Field* theField)
{
    assert(filename!=NULL);
    assert(theField!=NULL);

    FILE* theFile = fopen(filename, "w");
    fprintf(theFile, "<Field %d:Samples %d>\n", theField->fieldnum, theField->samples);
    int iPerm, jPerm;
    for(int i=0; i<theField->samples; i++)
    {
        for(int j=0; j<theField->samples; j++)
        {
            iPerm = theField->perm->index[i];
            jPerm = theField->perm->index[j];
            fprintf(theFile, "%f", theField->element[iPerm][jPerm]);
            if(j<theField->samples-1) fprintf(theFile, "\t");
        }
        fprintf(theFile, "\n");
    }
    fclose(theFile);
    
}
void saveFieldToHTML(const char* filename, Field* theField)
{
    assert(filename!=NULL);
    assert(theField!=NULL);

    FILE* theFile = fopen(filename, "w");
    fprintf(theFile, "<table>\n");
    int iPerm, jPerm;
    for(int i=0; i<theField->samples; i++)
    {
        fprintf(theFile, "<tr>");
        for(int j=0; j<theField->samples; j++)
        {
            iPerm = theField->perm->index[i];
            jPerm = theField->perm->index[j];
            fprintf(theFile, "<td>%f</td>", theField->element[iPerm][jPerm]);
        }
        fprintf(theFile, "</tr>\n");
    }
    fprintf(theFile, "</table>");
    fclose(theFile);
}

void displayField(Field* theField)
{
    assert(theField!=NULL);
    int index=0;
    float currentElement;
    int iPerm, jPerm;
    for(int x=0; x<theField->samples; x++)
    {
        for(int y=0; y<theField->samples; y++)
        {
            iPerm = theField->perm->index[x];
            jPerm = theField->perm->index[y];
            if(x!=y)
            {
                currentElement = theField->element[iPerm][jPerm];
                printf("(%d,%d)%3.2f{%2d} ", iPerm, jPerm, currentElement, index);
                index++;
            } else {
                printf("(%d,%d)  ><     ", iPerm, jPerm);
            }
        }
        printf("\n");
    }
}


#pragma mark Landscapes

Landscape* allocateLandscape(void)
{
    return malloc(sizeof(Landscape));
}

Landscape* makeLandscapeFromTDVs(int files, const char* filenames[])
{
    assert(files>0);
    assert(filenames!=NULL);
    
    Landscape* theScape = NULL;
    theScape = allocateLandscape();
    theScape->numFields = files;
    theScape->fields = allocateArrayOfFields(files);
    theScape->numNonDiagElts=0;
    int n;
    
    for(int i=0; i<files; i++)
    {
        theScape->fields[i] = makeFieldFromTDV(filenames[i]);
        n=theScape->fields[i]->samples;
        theScape->numNonDiagElts  += n*n-n;
    }
    theScape->isRaw = true;
    theScape->isRanked=false;
    theScape->isRankBased=false;
    theScape->isCentered=false;
    theScape->hasFlatVersion=false;
    theScape->flatVersion=NULL;
    
    return theScape;
}

bool isInLowerDiagonal(int i, int j);
bool isInLowerDiagonal(int i, int j){ return i<j; }

void modifyLandscapeMeanify(Landscape* theData)
{
    assert(theData!=NULL);
    
    float theTotal = 0.0;
    int theCount=0;
    
    //Has this already been done?
    if(theData->isCentered) return;
    
    //Processed data isn't raw
    theData->isRaw=false;

    //Is this formerly-ranked data?
    if(theData->isRanked)
    {
        theData->isRankBased=true;
    } else {
        theData->isRankBased=false;
    }

    //Either way, it's not ranked anymore.
    theData->isRanked=false;
    
    
    //old flat version now invalid. Clear out the memory it used.
    theData->hasFlatVersion=false; 
    if(theData->flatVersion != NULL)
    {
        free(theData->flatVersion);
        theData->flatVersion=NULL;
    }
    
    int currentSampleSize;
    
    //Step one: Find mean
    for(int f=0; f<theData->numFields; f++)
    {
        currentSampleSize = theData->fields[f]->samples;
        
        for(int x=0; x<currentSampleSize; x++)
        {
            for(int y=0; y<currentSampleSize; y++)
            {
                //Caveat: Skip main diagonal
                if(x!=y) 
                {
                    //2x efficiency: Only need upper triangle
                    if(!(UPPER_ONLY && isInLowerDiagonal(x,y)))
                    {
                        theTotal += theData->fields[f]->element[x][y];
                        theCount++;
                    }
                }
            }
        }
    }
    float theMean = theTotal / theCount;
    
    //Step two: subtract off the mean
    for(int f=0; f<theData->numFields; f++)
    {
        currentSampleSize = theData->fields[f]->samples;
        
        for(int x=0; x<currentSampleSize; x++)
        {
            for(int y=0; y<currentSampleSize; y++)
            {
                //Caveat: Skip main diagonal and one redundant triangle
                if(x!=y && x>y) 
                {
                    theData->fields[f]->element[x][y] -= theMean;
                    theData->fields[f]->element[y][x] -= theMean;
                }
            }
        }
    }
    //It's now centered!
    theData->isCentered=true;
}

//CHANGE to operation on Landscapes (buncha fields)
void  modifyLandscapeRankify(Landscape* theData)
{
    assert(theData!=NULL);
    
    //Make sure this isn't redundant effort
    if(theData->isRanked) return;
    
    //Processed data isn't raw
    theData->isRaw=false;
    
    //Ranking uncenters data
    theData->isCentered=false;
    

    List* sorted = makeListFromLandscape(theData); 
    modifyListSortify(sorted);
    //Old flat version is unranked and out of order, so nix it
    theData->hasFlatVersion=false; 
    theData->flatVersion=NULL;
    
    rankAndCount* rankInfo = NULL;
    
    int currentSampleSize;
    
    for(int f=0; f<theData->numFields; f++)
    {
        currentSampleSize = theData->fields[f]->samples;
        
        for(int x=0; x<currentSampleSize; x++)
        {
            for(int y=0; y<currentSampleSize; y++)
            {
                //Caveat: Skip main diagonal and one (redundant) triangle
                if(x!=y && x>y) 
                {
                    rankInfo = computeRankInList(theData->fields[f]->element[x][y],
                                                 sorted, rankInfo);
                    theData->fields[f]->element[x][y] = rankInfo->rank;
                    theData->fields[f]->element[y][x] = rankInfo->rank;
                }
            }
        }
    }

    free(rankInfo);
    //sorted list should no longer be in use by anyone
    free(sorted);
    
    //We're now ranked!
    theData->isRanked=true;
    //Ranked data is also rank-based
    theData->isRankBased=true;
}


#pragma mark Correlation
CorrelationAggregate* allocateCA(void)
{
    return malloc(sizeof(CorrelationAggregate));
}

void initializeCA(CorrelationAggregate* theCA)
{
    assert(theCA!=NULL);
    
    theCA->numerator = 0;
    theCA->denominatorL = 0;
    theCA->denominatorR = 0;
}
void augmentCAByValues(CorrelationAggregate* theCA, float xmxbar, float ymybar)
{
    assert(theCA!=NULL);
    
    theCA->numerator += xmxbar*ymybar;
    theCA->denominatorL += xmxbar*xmxbar;
    theCA->denominatorR += ymybar*ymybar;
}
float finishCorrelation(CorrelationAggregate* theCA)
{
    assert(theCA!=NULL);
    
    float n,l,r;
    n=theCA->numerator;
    l=theCA->denominatorL;
    r=theCA->denominatorR;
    return (FLOATIFY*n)/(sqrt(l*r));
}

void augmentCAByFields(CorrelationAggregate* theCA, Field* X, Field* Y)
{
    assert(X!=NULL);
    assert(Y!=NULL);
    assert(theCA!=NULL);
    assert(X->samples == Y->samples);
    
    //Compute correlation
    
    //initializeCA(theCA);
    float yVal;
    float xVal;
    int iPerm, jPerm;
    
    
    if(VERBOSE)
    {
        printf("Aggregating fields: \n");
        displayField(X);
        printf("AND\n");
        displayField(Y);
        printf("CA: N=%f, L=%f, R=%f\n", theCA->numerator,
               theCA->denominatorL, theCA->denominatorR);
    }
    
    for(int i=0; i<X->samples; i++)
    {
        if(VERBOSE) printf("\ti=%d\n",i);
        for(int j=0; j<X->samples; j++)
        {
            if(VERBOSE) printf("\t\tj=%d",j);
            //Caveat: Skip main diagonal 
            if(i!=j)
            {
                iPerm = Y->perm->index[i];
                jPerm = Y->perm->index[j];
                yVal  = Y->element[iPerm][jPerm];

                iPerm = X->perm->index[i];
                jPerm = X->perm->index[j];
                xVal  = X->element[iPerm][jPerm];
                //2x efficiency: Only need upper triangle
                if(!(UPPER_ONLY && isInLowerDiagonal(i, j)))
                {
                    if(VERBOSE) printf("\t\t\tx=%f, y=%f, ", xVal, yVal);
                    augmentCAByValues(theCA, xVal, yVal);
                    if(VERBOSE) printf("n=%f, L=%f, R=%f", theCA->numerator, 
                           theCA->denominatorL, theCA->denominatorR);
                }
            }
            if(VERBOSE) printf("\n");
        }
    }
}

float mantelR(Landscape* mPreserved, 
              Landscape* mPermuted,
              CorrelationAggregate* theCA)
{
    assert(mPermuted!=NULL);
    assert(mPreserved!=NULL);
    assert(mPermuted->numFields==mPreserved->numFields);
    assert(mPreserved->isCentered);
    assert(mPermuted->isCentered);
    
    //If no aggregator supplied, create a temporary one.
    bool noAggregator = (theCA == NULL);
    if(noAggregator) theCA = allocateCA();
    
    //Compute partial sums
    initializeCA(theCA);
    for(int f=0; f<mPermuted->numFields; f++)
    {
        augmentCAByFields(theCA, mPreserved->fields[f], mPermuted->fields[f]);
    }
    
    //Process aggregated correlator
    float theCorrelation = finishCorrelation(theCA);
    if(noAggregator) free(theCA);
    return theCorrelation;
}

float mantelRPartial(Landscape* mPermuted,
                        Landscape* mPreserved, 
                        Landscape* mGiven,
                        CorrelationAggregate* theCA)
{
    //If no aggregator supplied, create a temporary one.
    bool noAggregator = (theCA == NULL);
    if(noAggregator) theCA = allocateCA();

    float AB = mantelR(mPreserved, mPermuted, theCA);
    float AC = mantelR(mGiven,     mPermuted, theCA);
    float BC = mantelR(mPreserved, mGiven,    theCA);
    
    if(noAggregator) free(theCA);

    return (AB-AC*BC)/sqrt((1-AC*AC)*(1-BC*BC));    
}

StatisticalData* allocateStatData(void)
{
    return malloc(sizeof(StatisticalData));
}

StatisticalData* correlateAndFindP(Landscape* lPermuted, 
                                   Landscape* lPreserved, 
                                   int trials)
{
    assert(lPermuted!=NULL);
    assert(lPreserved!=NULL);
    assert(lPermuted->numFields == lPreserved->numFields);
    StatisticalData* theResults=allocateStatData();
    theResults->listOfCorrelations = allocateList();
    theResults->listOfCorrelations->count = trials;
    theResults->listOfCorrelations->data = allocateArrayOfFloats(trials);
    float currentCor;
    CorrelationAggregate* aCA = allocateCA();

    theResults->correlationType = "Unset";

    modifyLandscapeMeanify(lPreserved);
    modifyLandscapeMeanify(lPermuted);

    for(int perm=0; perm<trials; perm++)
    {
        //Make new random permutations
        for(int f=0; f<lPreserved->numFields; f++)
        {
            //Identity permutation the first time through, random the rest
            modifyPermPermutify(lPermuted->fields[f]->perm,
                                (!perm)?SEED_IDENTITY:SEED_RANDOM);
        }
        
        currentCor=mantelR(lPreserved, lPermuted, aCA);
        theResults->listOfCorrelations->data[perm] = currentCor;
        //Store the first result specially
        if(!perm) theResults->correlationOfInterest = currentCor;
    }
    modifyListSortify(theResults->listOfCorrelations);
    theResults->rankInfo = computeRankInList(theResults->correlationOfInterest, 
                                             theResults->listOfCorrelations, 
                                             NULL);
    return theResults;
}

StatisticalData* correlatePartialAndFindP(Landscape* lPermuted, 
                                          Landscape* lPreserved, 
                                          Landscape* lGiven, 
                                          int trials)
{
    assert(lPermuted!=NULL);
    assert(lPreserved!=NULL);
    assert(lGiven!=NULL);
    
    assert(lPermuted->numFields == lPreserved->numFields);
    assert(lGiven->numFields == lPreserved->numFields);
    
    StatisticalData* theResults=allocateStatData();
    theResults->listOfCorrelations = allocateList();
    theResults->listOfCorrelations->count = trials;
    theResults->listOfCorrelations->data = allocateArrayOfFloats(trials);
    float currentCor;
    CorrelationAggregate* aCA = allocateCA();
    
    theResults->correlationType = "Unset";
    
    modifyLandscapeMeanify(lPreserved);
    modifyLandscapeMeanify(lPermuted);
    modifyLandscapeMeanify(lGiven);
    
    for(int perm=0; perm<trials; perm++)
    {
        //Make new random permutations
        for(int f=0; f<lPreserved->numFields; f++)
        {
            //Identity permutation the first time through, random the rest
            modifyPermPermutify(lPermuted->fields[f]->perm,
                                (!perm)?SEED_IDENTITY:SEED_RANDOM);
        }
        
        currentCor=mantelRPartial(lPermuted, lPreserved, lGiven, aCA);
        theResults->listOfCorrelations->data[perm] = currentCor;
        //Store the first result specially
        if(!perm) theResults->correlationOfInterest = currentCor;
    }
    modifyListSortify(theResults->listOfCorrelations);
    theResults->rankInfo = computeRankInList(theResults->correlationOfInterest, 
                                             theResults->listOfCorrelations, 
                                             NULL);
    return theResults;
}


void saveData(StatisticalData* dataToSave, int timestamp)
{
    assert(dataToSave!=NULL);
    
    FILE* output = NULL;
    char fname[100];
    
    int trials = dataToSave->listOfCorrelations->count;
    float rank = dataToSave->rankInfo->rank+1;
    
    sprintf(fname, "testinfo.%d.%s.%s", 
            timestamp, dataToSave->correlationType, "txt");
    output = fopen(fname, "w");
    fprintf(output, "%s correlation is %f\n", 
            dataToSave->correlationType,
            dataToSave->correlationOfInterest);
    fprintf(output, "In %d trials, %d matches had average rank %f\n", 
            trials, dataToSave->rankInfo->count, rank);
    fprintf(output, "Assuming uniqueness, %f of the trials were >=\n", 
            (trials - rank)/(FLOATIFY*trials));
    fprintf(output, "Assuming uniqueness, %f of the trials were <=\n", 
            (rank)/(FLOATIFY*trials));
    fclose(output);
    
    sprintf(fname, "testinfo.%d.%s.%s", 
            timestamp, dataToSave->correlationType, "tdv");
    saveListToTDV(fname, dataToSave->listOfCorrelations);
}

void processFilePairs(int trials, int filesets, const char* argv[], int timestamp)
{
    const char* s[filesets]; //static: distances
    const char* p[filesets]; //permuted: differences
    int groupSize=2;
    int firstFile=2;
    for(int i=0; i<filesets; i++)
    {
        s[i] = argv[groupSize*i+firstFile];
        p[i] = argv[groupSize*i+firstFile+1];
    }
    
    //Load data
    Landscape* lPreserved = makeLandscapeFromTDVs(filesets, s);
    Landscape* lPermuted = makeLandscapeFromTDVs(filesets, p);
    
    StatisticalData* theStats = NULL;
    
    //Pearson correlation
    theStats = correlateAndFindP(lPermuted, lPreserved, trials);
    theStats->correlationType = "Pearson";
    saveData(theStats, timestamp);
    
    //Rank data
    modifyLandscapeRankify(lPreserved);
    modifyLandscapeRankify(lPermuted);
    
    //Spearman correlation
    theStats = correlateAndFindP(lPermuted, lPreserved, trials);
    theStats->correlationType = "Spearman";
    saveData(theStats, timestamp);
}

void processFileTriples(int trials, int filesets, const char* argv[], int timestamp)
{
    const char* s[filesets]; //static: distances
    const char* p[filesets]; //permuted: differences
    const char* g[filesets]; //given: differences to relativize by
    int groupSize = 3;
    int firstFile = 2;
    for(int i=0; i<filesets; i++)
    {
        s[i] = argv[groupSize*i+firstFile];
        p[i] = argv[groupSize*i+firstFile+1];
        g[i] = argv[groupSize*i+firstFile+2];
    }
    
    //Load data
    Landscape* lPreserved = makeLandscapeFromTDVs(filesets, s);
    Landscape* lPermuted = makeLandscapeFromTDVs(filesets, p);
    Landscape* lGiven = makeLandscapeFromTDVs(filesets, g);
    
    StatisticalData* theStats = NULL;
    
    //Pearson correlation
    theStats = correlatePartialAndFindP(lPermuted, lPreserved, lGiven, trials);
    theStats->correlationType = "Pearson (Partial)";
    saveData(theStats, timestamp);
    
    //Rank data
    modifyLandscapeRankify(lPreserved);
    modifyLandscapeRankify(lPermuted);
    modifyLandscapeRankify(lGiven);
    
    //Spearman correlation
    theStats = correlatePartialAndFindP(lPermuted, lPreserved, lGiven, trials);
    theStats->correlationType = "Spearman (Partial)";
    saveData(theStats, timestamp);
}


#pragma mark Lists
List* allocateList(void)
{
    return malloc(sizeof(List));
}

List* makeListFromList(List* theData)
{
    assert(theData!=NULL);
    List* theCopy = allocateList();
    theCopy->isSorted = theData->isSorted;
    theCopy->mean = theData->mean;
    theCopy->isMeanValid = theData->isMeanValid;
    theCopy->count = theData->count;
    theCopy->data = allocateArrayOfFloats(theCopy->count);
    for(int i=0; i<theCopy->count; i++)
    {
        theCopy->data[i] = theData->data[i];
    }
    return theCopy;
}

void  displayList(List* theData)
{     
    assert(theData!=NULL);
    printf("List %p: [ ", theData);
    for(int i=0; i<theData->count; i++)
    {
        printf("%.*f ", DISPLAY_DIGITS, theData->data[i]);
    }
    printf("]\n");
}

void  saveListToTDV(const char* filename, List* theData)
{
    assert(theData!=NULL);
    FILE* fout = fopen(filename, "w");
    assert(fout!=NULL);
    for(int i=0; i<theData->count; i++)
    {
        fprintf(fout, "%f\t", theData->data[i]);
    }
    fclose(fout);
}

void modifyListMeanify(List* theData)
{
    assert(theData!=NULL);
    
    float theTotal = 0.0;
    int theCount=0;
    
    if(!theData->isMeanValid)
    {
        for(int i=0; i<theData->count; i++)
        {
            theTotal += theData->data[i];
            theCount++;
        }
        theData->mean = theTotal / theCount;
        theData->isMeanValid = true;
    }
    
    for(int i=0; i<theData->count; i++)
    {
        theData->data[i] -= theData->mean;
    }

    theData->mean = 0.0;
}

////////////////////////////////////////////////////
// Heapsort adapted from 
// http://www.algorithmist.com/index.php/Heap_sort.c

void  modifyListSortify(List* theData)
{
    assert(theData!=NULL);
    
    if(theData->isSorted) return; //Why bother? It's sorted already
    
    for (int i = (theData->count / 2); i >= 0; i--) 
        heapSortSiftDown(theData, i, theData->count-1);
    
    for (int i = theData->count-1; i >= 1; i--) 
    {
        swapF(&theData->data[0], &theData->data[i]);
        heapSortSiftDown(theData, 0, i-1);
    }
    
    theData->isSorted = true; //Flag result for future timesaving
}

void  heapSortSiftDown(List* theData, int root, int bottom)
{
    assert(theData!=NULL);
    
    //This process disorders things
    theData->isSorted = false;

    int farChild = root * 2 + 1;
    int sibling;
    
    while(farChild <= bottom)
    {
        // Find the biggest child
        if(farChild < bottom) 
        {
            sibling = farChild + 1;
            // Reversed for stability
            if(!inOrder(theData->data[farChild], theData->data[sibling]))
            {
                farChild = sibling;
            }
        }
        
        // If we have the correct ordering, we are done.
        if(inOrder(theData->data[root], theData->data[farChild])) return;
        
        // Swap
        swapF(&theData->data[root], &theData->data[farChild]);
        
        root=farChild;
        farChild = root * 2 + 1;
    }
}

rankAndCount* allocateRankAndCount(void)
{
    return malloc(sizeof(rankAndCount));
}
rankAndCount* computeRankInList(float datum, List* theData,
                                        rankAndCount* theRank)
{
    assert(theData!=NULL);
    
    bool noRank = (theRank == NULL);
    if(noRank) theRank = allocateRankAndCount();
    
    bool workingOnCopy = ! theData->isSorted;
    
    if(workingOnCopy) 
    {
        printf("Note: Your data is being copied to sort.\n");
        theData = makeListFromList(theData);
        modifyListSortify(theData);   
    }
    
    //Binary search for a match
    int iguess = (theData->count)/2; //Search index
    int ilo=0;  //Low index
    int ihi=theData->count; //High index
    float candidate = theData->data[iguess]; //Found data
    
    while(candidate != datum) 
    {
        //DELETEME printf("Datum(%f)!=Candidate(%f). <%d | %d | %d>\n", datum, candidate, ilo, iguess, ihi);
        if(ihi<ilo)
        { //Candidate does not exist in list
            theRank->count = 0; 
            theRank->rank = (ihi+ilo)/2.0;
            return theRank;
        }
        //DELETEME printf("DBEUG: Check this faster search work!\n");
        if(inOrder(candidate, datum)) ihi=iguess-1;
        else ilo=iguess+1;
        iguess = (ilo+ihi)/2;
        candidate = theData->data[iguess];
    }
    
    ilo=ihi=iguess; //Expand match to range
    while(theData->data[ilo]==datum && ilo>=0) ilo--; //find pre-firt-match
    while(theData->data[ihi]==datum && ihi<theData->count) ihi++; //find post-last-match
    
    //Construct answer
    theRank->count = ((ihi-1)-ilo); //Elements between low and high
    theRank->rank = (ilo+ihi)/2.0; //Average position btwn low and high
    theRank->indexEnd = ihi-1;
    theRank->indexStart = ilo+1;
    
    if(workingOnCopy){ free(theData); } //No memory leak!
    
    return theRank;    
}

#pragma mark ->List

List*  makeListFromLandscape(Landscape* theData)
{
    assert(theData!=NULL);
    if(theData->hasFlatVersion) return theData->flatVersion;
    
    //Create a list of the appropriate size
    List* theList = allocateList();
    theList->count = theData->numNonDiagElts;
    theList->data = allocateArrayOfFloats(theList->count);
    theList->isSorted = false;
    
    //Fill list with data
    float currentElement;
    Field* theField;
    int index=0;
    float runningTotal=0.0;

    for(int f=0; f<theData->numFields; f++)
    {
        theField = theData->fields[f];
        
        for(int x=0; x<theField->samples; x++)
        {
            for(int y=0; y<theField->samples; y++)
            {
                //Caveat: Skip main diagonal
                if(x!=y) 
                {
                    index++;
                    currentElement = theField->element[x][y];
                    theList->data[index] = currentElement;
                    runningTotal += currentElement;
                }
            }
        }
    }
    //May as well have the mean on hand, just in case.
    theList->mean = runningTotal / index;
    theList->isMeanValid = true;
    
    //Store list for later re-use
    theData->hasFlatVersion = true;
    theData->flatVersion = theList;
    
    return theList;
}

