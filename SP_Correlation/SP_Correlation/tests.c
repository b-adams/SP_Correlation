/**
 * @file tests.c
 * @author Bryant Adams
 * @date 3/20/12
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <stdio.h>
#include "tests.h"
//#include "defines.h"
#include "functions.h"

#define TEST_SEED 31415
void makeTestFiles(void)
{
    Field* tmp;
    char filename[20];
    int size;
    for(int fieldNum=1; fieldNum<4; fieldNum++)
    {
        size = 1+fieldNum;
        for(int a=1; a<=2; a++)
        {
            tmp = makeRandomField(size);
            sprintf(filename, "Test%d%c.tdv", fieldNum, (a==1?'L':'R'));
            saveFieldToTDV(filename, tmp);
            free(tmp);
            tmp=NULL;
        }
    }
}

bool runTests(void)
{
    printf("Runtests still stubbed\n");
    void makeTestFiles();
    
    assert(testRandInRange());

    assert(testInOrder());

    assert(testSwapF());

    assert(testSwapI());

    assert(testInitializeCA());
 
    assert(testAugmentCAByValues());

    assert(testFinishCorrelation());

    assert(testMakePerm());

    assert(testModifyPermPermutify());

    assert(testMakeListFromList());

    //    assert( testDisplayList());

//    assert( testSaveListToTDV());

    assert( testModifyListMeanify());

    assert( testModifyListSortify());

    assert(testComputeRankInList());

    assert(testMakeRandomField());

    assert(testMakeFieldFromTDV());

    assert( testSaveFieldToTDV());

//    assert(testDisplayField());

#warning tests unimplemented
    assert(testAugmentCAByFields());
   
    assert(testMakeLandscapeFromTDVs());
    
    assert(testModifyLandscapeMeanify());
    
    assert(testModifyLandscapeRankify());
   
    assert(testMakeListFromField());
    
    assert(testMakeListFromLandscape());
    
    assert(testMantelR());
    assert(testMantelRPartial());
    
    assert(testCorrelateAndFindP());
    
    assert(testCorrelatePartialAndFindP());
    
    assert(testProcessFilePairs());
    
    assert(testProcessFileTriples());
    
    return reportEnd(true, NULL);
    
}

void reportStart(char* function)
{
    printf("Testing %s: ", function);
}
bool reportEnd(bool good, char* explanation)
{
    if(good)
    {
        printf("OK\n");
    } else {
        printf("BAD(%s)\n", explanation);
    }
    return good;
}
bool reportFake(void)
{
    printf("Temporarily acceptible\n");
    return true;
}
bool testRandInRange(void)
{
    reportStart("randInRange");
    int numtests=10000;
    int rand;
    int low, high;
    low = 14;
    high = 99;
    for(int i=0; i<numtests; i++)
    {
        rand = randInRange(low, high);
        if(rand<low) return reportEnd(false, "too low");
        if(rand>high) return reportEnd(false, "too high");
    }
    return reportEnd(true, NULL);
}

bool testInOrder(void)
{
    reportStart("inOrder");
    float before=5.0; //Currently sorting high-to-low
    float after=3.0;
    if(!inOrder(before, after)) return reportEnd(false, "false negative");
    if(inOrder(after, before)) return reportEnd(false, "false positive");
    return reportEnd(true, NULL);
}

bool testSwapF(void)
{
    reportStart("swapF");
    float lbefore,lafter,rbefore,rafter;
    lbefore = 3.0;
    rbefore = 5.0;
    lafter = lbefore;
    rafter = rbefore;
    swapF(&lafter, &rafter);
    if(lbefore != rafter) return reportEnd(false, "left not moved");
    if(rbefore != lafter) return reportEnd(false, "right not moved");
    return reportEnd(true, NULL);
}

bool testSwapI(void)
{
    reportStart("swapI");
    int lbefore,lafter,rbefore,rafter;
    lbefore = 3;
    rbefore = 5;
    lafter = lbefore;
    rafter = rbefore;
    swapI(&lafter, &rafter);
    if(lbefore != rafter) return reportEnd(false, "left not moved");
    if(rbefore != lafter) return reportEnd(false, "right not moved");
    return reportEnd(true, NULL);
}

#pragma mark Correlation

bool testInitializeCA(void)
{
    reportStart("initializeCA");
    CorrelationAggregate* pCA = allocateCA();
    initializeCA(pCA);
    if(pCA->numerator != 0) return reportEnd(false, "nonzero numerator");
    if(pCA->denominatorL != 0) return reportEnd(false, "nonzero Ldenom");
    if(pCA->denominatorR != 0) return reportEnd(false, "nonzero Rdenom");
    return reportEnd(true, NULL);
}

bool testAugmentCAByValues(void)
{
    assert(testInitializeCA());
    reportStart("augmentCAByValues");
    CorrelationAggregate* theCA = allocateCA();
    initializeCA(theCA);
    float xmxb = 3.0;
    float ymyb = 5.0;
    
    augmentCAByValues(theCA, xmxb, ymyb);
    if(theCA->numerator != xmxb*ymyb) return reportEnd(false, "numerator miscalc");
    if(theCA->denominatorL != xmxb*xmxb) return reportEnd(false, "X denominator miscalc");
    if(theCA->denominatorR != ymyb*ymyb) return reportEnd(false, "Y denominator miscalc");
        
    return reportEnd(true, NULL);
}


bool testFinishCorrelation(void)
{
    assert(testInitializeCA());
    assert(testAugmentCAByValues());
    reportStart("finishCorrelation");
    
    CorrelationAggregate* theCA = allocateCA();
    initializeCA(theCA);
    float x1=3.0;
    float y1=4.0;
    float x2=4.0;
    float y2=3.0;
    augmentCAByValues(theCA, x1,y1);
    augmentCAByValues(theCA, x2,y2);
    float correlationValue = finishCorrelation(theCA);
    float expectedResult = (x1*y1+x2*y2)/sqrt((x1*x1+x2*x2)*(y1*y1+y2*y2));
    if(fabs(correlationValue - expectedResult) > 0.00001)
        return reportEnd(false, "incorrect result");  
    return reportEnd(true, NULL);
}


#pragma mark Permutation

bool testMakePerm(void)
{
    reportStart("makePerm");
    int testSize=3;
    Perm* aPerm = makePerm(testSize, TEST_SEED);
    if(aPerm->size != testSize) return reportEnd(false, "wrong size");
    int* indices_found = allocateArrayOfInts(testSize);
    for(int i=0; i<testSize; i++)
    {
        indices_found[i]=0;
    }
    for(int i=0; i<testSize; i++)
    {
        indices_found[aPerm->index[i]]++;
    }
    for(int i=0; i<testSize; i++)
    {
        if(indices_found[i]<1)
            return reportEnd(false, "index skipped");
        if(indices_found[i]>1)
            return reportEnd(false, "index repeated");
    }
    
    aPerm = makePerm(testSize, SEED_IDENTITY);
    for(int i=0; i<testSize; i++)
    {
        if(aPerm->index[i] != i) return reportEnd(false, "identity fail");
    }

    return reportEnd(true, NULL);
}

void showPermutation(Perm* thePerm);
void showPermutation(Perm* thePerm)
{
    printf("<Perm>[ ");
    for(int i=0; i<thePerm->size; i++)
    {
        printf("%d->%d ",i, thePerm->index[i]);
    }
    printf("]\n");
}
bool testModifyPermPermutify(void)
{
    assert(testMakePerm());
    reportStart("modifyPermPermutify");
    int testSize = 5;
    Perm* thePerm = makePerm(testSize, SEED_IDENTITY);
    //showPermutation(thePerm);
    for(int i=0; i<testSize; i++)
    {
        if(thePerm->index[i] != i) return reportEnd(false, "pre-identity fail");
    }
    for(int i=0; i<testSize-1; i++)
    {
        modifyPermPermutify(thePerm, TEST_SEED);
        //showPermutation(thePerm);
    }   //Should result in identity
    for(int i=0; i<testSize; i++)
    {
        if(thePerm->index[i] != i) return reportEnd(false, "post-identity fail");
    }
    modifyPermPermutify(thePerm, TEST_SEED); //Mess it up once
    int displacements=0;
    for(int i=0; i<testSize; i++)
    {
        if(thePerm->index[i] != i) displacements++;
    }
    if(displacements<1) return reportEnd(false, "nonidentity fail");
    
    return reportEnd(true, NULL);
}


#pragma mark Lists

List* makeTestList(void);
List* makeTestList(void)
{
    List* theList = allocateList();
    int testEntries = 12;
    theList->count = testEntries;
    theList->data = allocateArrayOfFloats(testEntries);
    for(int i=0; i<testEntries; i++)
    {
        theList->data[i] = (randInRange(1, 99)/10.0);
    }
    theList->isMeanValid = false;
    theList->isSorted = false;
    return theList;
}

bool testMakeListFromList(void)
{
    reportStart("makeListFromList");

    List* oldData = makeTestList();
    
    List* newData = makeListFromList(oldData);
    //displayList(oldData);
    //displayList(newData);
    
    if(oldData->count != newData->count) return reportEnd(false, "size mismatch");
    for(int i=0; i<oldData->count; i++)
    {

        if(fabs(oldData->data[i] - newData->data[i]) > 0.0001)
            return reportEnd(false, "data mismatch");
    }
    if(oldData->isSorted != newData->isSorted) return reportEnd(false, "sortV mismatch");
    if(oldData->isMeanValid != newData->isMeanValid) return reportEnd(false, "meanV mismatch");
    if(oldData->isMeanValid && (oldData->mean != newData->mean)) return reportEnd(false, "mean mismatch");
          
    return reportEnd(true, NULL);
}


bool  testDisplayList(void)
{
    reportStart("displayList");
    displayList(makeTestList());
    printf("USER: did you just see a list?");
    char answer='w';
    while(answer!='Y' && answer!='N')
    {
        printf(" Y/N: ");
        scanf(" %c", &answer);
    }
    return reportEnd(answer=='Y', "User rejected");
}


bool  testSaveListToTDV(void)
{
    reportStart("saveListToTDV");
    saveListToTDV("testSave.tdv", makeTestList());
    
    printf("USER: please check testSave.tdv. Is it ok?");
    char answer='w';
    while(answer!='Y' && answer!='N')
    {
        printf(" Y/N: ");
        scanf(" %c", &answer);
    }
    return reportEnd(answer=='Y', "User rejected");
}


bool  testModifyListMeanify(void)
{
    reportStart("modifyListMeanify");
    List* theData = makeTestList();
    modifyListMeanify(theData);
    float theSum = 0;
    if(!(theData->isMeanValid)) return reportEnd(false, "mean not valid");
    if(fabs(theData->mean) > 0.0001) return reportEnd(false, "mean not zero");
    
    for(int i=0; i<theData->count; i++)
    {
        theSum += theData->data[i];
    }
    
    if(fabs(theSum) > 0.0001) return reportEnd(false, "mean miscalculated");
    return reportEnd(true, NULL);
}


bool  testModifyListSortify(void)
{    
    reportStart("modifyListSortify");
    List* theData = makeTestList();
//    displayList(theData);
    modifyListSortify(theData);
//    displayList(theData);
    if(!(theData->isSorted)) return reportEnd(false, "sort not valid");
    if(fabs(theData->mean) > 0.0001) return reportEnd(false, "mean not zero");
    float lo, hi;
    
    for(int i=0; i<theData->count-1; i++)
    {
        lo = theData->data[i];
        hi = theData->data[i+1];
        //Sorting algorithm gives lists in reverse order
        if(!inOrder(hi, lo))
        {
            return reportEnd(false, "out of order");
        }
    }
    return reportEnd(true, NULL);
}


bool testComputeRankInList(void)
{   reportStart("computeRankInList");
    float datum;
    List* theData;
    rankAndCount* theRank;
    theData = makeTestList();
    datum = theData->data[theData->count / 2];
    theRank = computeRankInList(datum, theData, NULL);
    
    modifyListSortify(theData);
//    displayList(theData);
    if(theData->data[(int)theRank->rank] != datum)
        return reportEnd(false, "incorrect rank");
    int countCheck = 0;
    for(int i=0; i<theData->count; i++)
    {
        if(theData->data[i] == datum) 
        {
            
        countCheck++;
            if(i < theRank->indexStart) 
                reportEnd(false, "bad start");
            if(i > theRank->indexEnd) 
                reportEnd(false, "bad end");
        } else {
            if(theRank->indexStart < i && i < theRank->indexEnd) 
                reportEnd(false, "nondatum in range");
        }

    }
    if(countCheck!=theRank->count) return reportEnd(false, "incorrect count");

    return reportEnd(true, NULL);
}



#pragma mark Fields

bool testMakeRandomField(void)
{
    reportStart("makeRandomField");
    int samples=10;
    Field* tfld = makeRandomField(samples);
    if(tfld->samples != samples) return reportEnd(false, "sample count");  
    for(int i=0; i<samples; i++)
    {
        for(int j=0; j<i; j++)
        {
            if(tfld->element[i][j] != tfld->element[j][i])
                return reportEnd(false, "not symmetric");
            
        }
    }
    return reportEnd(true, NULL);
}


bool testLoadAndSaveFieldTDV(void);
bool testLoadAndSaveFieldTDV(void)
{
    int samples = 10;
    Field* tfld = makeRandomField(samples);
    char* filename = "testMake.tdv";
    saveFieldToTDV(filename, tfld);
    Field* tfld2 = makeFieldFromTDV(filename);
    
    if(tfld->samples != tfld2->samples)
        return reportEnd(false, "sample size mismatch");
    
    if(tfld->fieldnum != tfld2->fieldnum)
        return reportEnd(false, "field num mismatch");
    
    for(int i=0; i<samples; i++)
    {
        for(int j=0; j<samples; j++)
        {
            if(tfld->element[i][j] != tfld2->element[i][j])
                return reportEnd(false, "element mismatch");
        }
    }

    return reportEnd(true, NULL);
}



bool testMakeFieldFromTDV(void)
{
    reportStart("makeFieldFromTDV");
    return testLoadAndSaveFieldTDV();
}
bool  testSaveFieldToTDV(void)
{
    reportStart("saveFieldToTDV");
    return testLoadAndSaveFieldTDV();
}


bool testDisplayField(void)
{
    reportStart("displayField");
    Field* td = makeRandomField(8);
    saveFieldToTDV("displayTest.txt", td);
    printf("\n");
    displayField(td);
    printf("USER: Does the displayed field match the contents of displayText.txt?");
    char answer='w';
    while(answer!='Y' && answer!='N')
    {
        printf(" Y/N: ");
        scanf(" %c", &answer);
    }
    return reportEnd(answer=='Y', "User rejected");
}




/**
 * @brief Augment a correlation aggregate with information from two fields
 * @param X First field to correlate
 * @param Y Second field to correlate
 * @param theCA Correlation aggregate to store cumulative information
 * @returns the correlation so far, which may not be the final result
 * @sideeffect Adds correlation information from fields X and Permute(Y) to theCA
 */
bool testAugmentCAByFields(void)
{//CorrelationAggregate* theCA, Field* X, Field* Y
    reportStart("augmentCAByFields");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}



#pragma mark Landscapes

bool testMakeLandscapeFromTDVs(void)
{//int files, const char* filename[]
    reportStart("makeLandscapeFromTDVs");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


bool testModifyLandscapeMeanify(void)
{//Landscape* theData
    reportStart("modifyLandscapeMeanify");
    return reportFake();
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


bool testModifyLandscapeRankify(void)
{//Landscape* theData
    reportStart("modifyLandscapeRankify");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


#pragma mark Field->List
/**
 * @brief Flatten a 2D field to a 1D list
 * @param theField Field data to flatten
 * @returns A List containing all the elements of theField
 * @sideeffect Adds flattened list to theField for memoization.
 */
bool testMakeListFromField(void)
{//Field* theField
    reportStart("makeListFromField");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


bool testMakeListFromLandscape(void)
{//Landscape* theData
    reportStart("makeListFromLandscape");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


#pragma mark Correlation
/**
 * @brief Find correlation of two sets of fields
 * @param theFieldSet DataSet containing base fields and permutation to use
 * @param theCA optional (NULL for auto) correlation aggregator
 * @returns correlation value between base fields
 */

bool testMantelR(void)
{//Landscape* mPreserved, Landscape* mPermuted, CorrelationAggregate* theCA
    reportStart("mantelR");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}

bool testMantelRPartial(void)
{//Landscape* mPermuted,    Landscape* mPreserved,     Landscape* mGiven,    CorrelationAggregate* theCA
    reportStart("mantelRPartial");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


#pragma mark P value

/**
 * @brief Run input data through multiple correlations
 * @param theFieldSet DataSet structure containing data to permute and correlate
 * @param trials Number of permutations to correlate
 * @returns StatisticalData on the rank of the first correlation among all the rest
 */
bool testCorrelateAndFindP(void)
{//Landscape* lPermuted,     Landscape* lPreserved,     int trials
    reportStart("correlateAndFindP");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


bool testCorrelatePartialAndFindP(void)
{//Landscape* lPermuted,     Landscape* lPreserved,     Landscape* lGiven,     int trials
    reportStart("correlatePartialAndFindP");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


/**
 * @brief Creates files containing data on Spearman and Person correlation of inputs
 * @param trials Number of permutations to correlate for each type
 * @param filepairs Number of substrata that will be supplied
 * @param argv Array of (original) command-line arguments.
 * @param timestamp Time used for random seed, and to put in filenames
 * @sideeffect Creates testinfo.TIMESTAMP.[Pearson|Spearman].[txt|tdv] files.
 */
bool testProcessFilePairs(void)
{//int trials, int filesets, const char* argv[], int timestamp
    reportStart("processFilePairs");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


bool testProcessFileTriples(void)
{//int trials, int filesets, const char* argv[], int timestamp
    reportStart("processFileTriples");
    return reportFake();
    return reportEnd(false, "STUB");  
    return reportEnd(true, NULL);
}


