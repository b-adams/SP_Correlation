/**
 * @file functions.h
 * @author Bryant Adams
 * @date 3/20/12
 */

#ifndef SC_functions_h
#define SC_functions_h

#include <string.h>
#include <stdbool.h>
#include "defines.h"

#pragma mark Utility

float* allocateArrayOfFloats(int howmany);
float** allocateArrayOfArraysOfFloats(int howmany);
int* allocateArrayOfInts(int howmany);

/**
 * @brief Generate random number in a given range
 * @param lo Smallest number allowable for output
 * @param hi Largest number allowable for output
 * @returns integer n, lo<=n<=hi
 */
int randInRange(int lo, int hi);

/**
 * @brief Check if inputs are in desired order
 * @param l Allegedly first element
 * @param r Allegedly second element
 * @returns TRUE when l<r, FALSE when r<l. 
 */
bool inOrder(float l, float r);

/**
 * @brief Exchange contents of two float variables
 * @param l A variable to swap
 * @param r A variable to swap
 * @sideeffect Exchanges contents of l and r
 */
void swapF(float* l, float* r);

/**
 * @brief Exchange contents of two int variables
 * @param l A variable to swap
 * @param r A variable to swap
 * @sideeffect Exchanges contents of l and r
 */
void swapI(int* l, int* r);

#pragma mark Correlation
/**
 * @brief Partial sum of components for a correlation computation
 */
typedef struct {
    float numerator;    /**< sum((x-xbar)*(y-ybar)) */
    float denominatorL; /**< sum((x-xbar)*(x-ybar)) */
    float denominatorR; /**< sum((y-xbar)*(y-ybar)) */
} CorrelationAggregate;

/**
 * @brief Initialize a correlation aggregate
 * @param theCA CorrelationAggregate to initialize
 * @sideeffect Sets all partial sums to zero
 */
void initializeCA(CorrelationAggregate* theCA);

CorrelationAggregate* allocateCA(void);

/**
 * @brief Add data to a correlation aggregate
 * @param theCA CorrelationAggregate to augment
 * @param xmxbar Value of x-xbar to use
 * @param ymybar Value of y-ybar to use
 * @sideeffect Adds (x-xbar)*(y-ybar) to numerator, similar quantities to denominators.
 */
void augmentCAByValues(CorrelationAggregate* theCA, float xmxbar, float ymybar);

/**
 * @brief Extract result from correlation aggregate
 * @param theCA CorrelationAggregate compute correlation from
 * @returns numerator/sqrt(denomL*denomR)
 */
float finishCorrelation(CorrelationAggregate* theCA);

#pragma mark Permutation
/**
 @brief Represents a permutation as a collection of index remappings
 */
typedef struct {
    int * index; /**< Array of indices */
    int size;    /**< Number of indices */
} Perm;

Perm* allocatePermutation(void);

/**
 * @brief Creates a permutation of a given size
 * @param size How many elements to permute
 * @param seed 0 for identity permutation, -1 for existing randomization, >0 for a known seed value
 * @returns a random permutation unless seed is 0
 */
Perm* makePerm(int size, int seed); //0 id, -1 no seed

/**
 * @brief Change a permutation
 * @param thePerm Permutation to permute
 * @param seed 0 for identity permutation, -1 for existing randomization, >0 for a known seed value
 * @sideeffect Permutes thePerm
 */
void modifyPermPermutify(Perm* thePerm, int seed);

#pragma mark Lists
/**
 * @brief Represents a vector of floats, tracking mean and sortedness.
 */
typedef struct {
    float * data;     /**< Array of floats */
    int count;        /**< Size of the array */
    float mean;       /**< Mean of all array entries */
    bool isMeanValid; /**< False unless mean has been calculated */
    bool isSorted;    /**< False unless list has been sorted */
} List;

List* allocateList(void);

/**
 * @brief Create a List by copying another list
 * @param theData List to copy
 * @returns a copy of theData
 */
List* makeListFromList(List* theData);

/**
 * @brief Display a list on screen
 * @param theData The list to display
 * @sideeffect Outputs contents of a list to screen
 */
void  displayList(List* theData);

/**
 * @brief Output a list to file
 * @param filename Name of file to store data in
 * @param theData The list to output
 * @sideeffect Create/overwrite 
 */
void  saveListToTDV(const char* filename, List* theData);

/**
 * @brief Center a list about its mean
 * @param theData The list to center
 * @sideeffect Set all elements X of theData to X-Xbar, where Xbar is mean of elements in theData 
 */
void  modifyListMeanify(List* theData);

/**
 * @brief Sort a list (heapsort)
 * @param theData The list to sort
 * @sideeffect theData will contain same elements, in sorted order 
 */
void  modifyListSortify(List* theData);

/**
 * @brief Utility function for sorting
 * @param theData The list bubble
 * @param root Node index to bubble down from
 * @param bottom Node index of last entry
 * @sideeffect theData will contain same elements, in sorted order 
 */
void  heapSortSiftDown(List* theData, int root, int bottom);

/**
 * @brief Represents positional information about an element in a list
 */
typedef struct {
    float rank; /**< Average position of all matching elements */
    int count;  /**< Number of elements matching */
    int indexStart;
    int indexEnd;
} rankAndCount;

rankAndCount* allocateRankAndCount(void);

rankAndCount* computeRankInList(float datum, List* theData, rankAndCount* theRank);


#pragma mark Fields
/**
 * @brief Represents distance-matrix data from a field of samples.
 */
typedef struct {
    int fieldnum;        /**< User name for field. Not really used. */
    int samples;         /**< Number of samples being compared */
    float ** element;    /**< 2D array (samples^2 elements) of comparisons */
//    bool isRaw;          /**< FALSE unless based on fresh data import. */
//    bool isRanked;       /**< FALSE unless this has been turned into ranked elements */
//    bool isRankBased;       /**< TRUE if ranked or ranked-and-mean'd */
//    bool isCentered;     /**< FALSE unless this has been centered about its mean */
    bool hasFlatVersion; /**< FALSE unless a flattened version has been generated */
    List* flatVersion;   /**< All elements arranged into a 1D array */
    Perm* perm;
} Field;

Field* allocateField(void);

/**
 * @brief Create a symmetric Field with random data
 * @param samples Width (and height) of the comparison matrix
 * @returns a randomly generated Field (symmetric, zero on diagonal)
 */
Field* makeRandomField(int samples);

/**
 * @brief Load a tab-delimited-value file into a Field
 * @param filename Filename for the TDV file.
 * @returns A Field with contents matching the TDV file.
 */
Field* makeFieldFromTDV(const char* filename);

/**
 * @brief Save field data to file.
 * @param filename Name of file to create/overwrite with data
 * @param theField Field structure to dump to file
 * @sideeffect creates/overwrites tab-delimited file
 */
void   saveFieldToTDV(const char* filename, Field* theField);

/**
 * @brief Save field data to file.
 * @param filename Name of file to create/overwrite with data
 * @param theField Field structure to dump to file
 * @sideeffect creates/overwrites an HTML Table file
 */
void   saveFieldToHTML(const char* filename, Field* theField);

/**
 * @brief Display field data on screen.
 * @param theField Field structure to dump to file
 * @sideeffect Displays field indices and flattened index around data
 */
void   displayField(Field* theField);

/**
 * @brief Augment a correlation aggregate with information from two fields
 * @param X First field to correlate
 * @param Y Second field to correlate
 * @param theCA Correlation aggregate to store cumulative information
 * @returns the correlation so far, which may not be the final result
 * @sideeffect Adds correlation information from fields X and Permute(Y) to theCA
 */
void augmentCAByFields(CorrelationAggregate* theCA, Field* X, Field* Y);


#pragma mark Landscapes
typedef struct {
    int numFields;
    int numNonDiagElts;
    Field** fields;
    bool isRaw;          /**< FALSE unless based on fresh data import. */
    bool isRanked;       /**< FALSE unless this has been turned into ranked elements */
    bool isRankBased;       /**< TRUE if ranked or ranked-and-mean'd */
    bool isCentered;     /**< FALSE unless this has been centered about its mean */
    bool hasFlatVersion; /**< FALSE unless a flattened version has been generated */
    List* flatVersion;   /**< All elements arranged into a 1D array */
} Landscape;

Landscape* allocateLandscape(void);

Landscape* makeLandscapeFromTDVs(int files, const char* filename[]);

void modifyLandscapeMeanify(Landscape* theData);

void  modifyLandscapeRankify(Landscape* theData);

#pragma mark Field->List
/**
 * @brief Flatten a 2D field to a 1D list
 * @param theField Field data to flatten
 * @returns A List containing all the elements of theField
 * @sideeffect Adds flattened list to theField for memoization.
 */
//List*  makeListFromField(Field* theField);

List*  makeListFromLandscape(Landscape* theData);

#pragma mark Correlation
/**
 * @brief Contains loaded/permuted data for analysis
 */
typedef struct {
    Landscape* lPreserved;
    Landscape* lPermuted;
    int numFields; /**< Number of fields being analyzed */
} DataSet;

/**
 * @brief Find correlation of two sets of fields
 * @param theFieldSet DataSet containing base fields and permutation to use
 * @param theCA optional (NULL for auto) correlation aggregator
 * @returns correlation value between base fields
 */

float mantelR(Landscape* mPreserved, 
                 Landscape* mPermuted,
                 CorrelationAggregate* theCA);
float mantelRPartial(Landscape* mPermuted,
                        Landscape* mPreserved, 
                        Landscape* mGiven,
                        CorrelationAggregate* theCA);

#pragma mark P value
/**
 * @brief Holds a value, a list of values, and information on that value's place in the list
 */
typedef struct {
    char* correlationType; /**< For when written to file (pearson or spearman) */
    float correlationOfInterest; /**< The value being held */
    rankAndCount* rankInfo; /**< Information about the value's place in the list */
    List* listOfCorrelations; /**< The list of other sample values */
} StatisticalData;

StatisticalData* allocateStatData(void);

/**
 * @brief Run input data through multiple correlations
 * @param theFieldSet DataSet structure containing data to permute and correlate
 * @param trials Number of permutations to correlate
 * @returns StatisticalData on the rank of the first correlation among all the rest
 */
StatisticalData* correlateAndFindP(Landscape* lPermuted, 
                                   Landscape* lPreserved, 
                                   int trials);

StatisticalData* correlatePartialAndFindP(Landscape* lPermuted, 
                                   Landscape* lPreserved, 
                                   Landscape* lGiven, 
                                   int trials);

/**
 * @brief Creates files containing data on Spearman and Person correlation of inputs
 * @param trials Number of permutations to correlate for each type
 * @param filepairs Number of substrata that will be supplied
 * @param argv Array of (original) command-line arguments.
 * @param timestamp Time used for random seed, and to put in filenames
 * @sideeffect Creates testinfo.TIMESTAMP.[Pearson|Spearman].[txt|tdv] files.
 */
void processFilePairs(int trials, int filesets, const char* argv[], int timestamp);

void processFileTriples(int trials, int filesets, const char* argv[], int timestamp);

/**
 * @brief Saves statistical data to file
 * @param timestamp Identifier to distinguish files from different runs
 * @param dataToSave Statistical data to output to file
 * @sideeffect Creates a .txt and a .tdv file containing summary info and data for a histogram.
 */
void saveData(StatisticalData* dataToSave, int timestamp);

#endif

//////// Useful links
// http://en.cppreference.com/w/cpp/io/c/fprintf
// http://www.gnu.org/savannah-checkouts/gnu/libc/manual/html_node/Freeing-after-Malloc.html
// http://www.cprogramming.com/tutorial/cfileio.html
// http://cboard.cprogramming.com/c-programming/131918-reading-tab-delimited-file.html
// http://www.cplusplus.com/reference/clibrary/cstring/strtok/
// http://www.cplusplus.com/reference/clibrary/cstdio/fgets/
// http://www.cplusplus.com/reference/clibrary/cstdio/sscanf/
// http://en.wikipedia.org/wiki/Fisher–Yates_shuffle
// http://en.wikipedia.org/wiki/Assert.h
// http://www.eggheadcafe.com/microsoft/VC-Language/33837208/const-char--vs-char-.aspx
// http://www.cplusplus.com/reference/clibrary/cstdlib/srand/
// http://stackoverflow.com/questions/6864308/newbie-random-number-generator-in-c-question
// http://www.cplusplus.com/reference/clibrary/cstdio/fscanf/
// http://www.cplusplus.com/reference/clibrary/cstdio/sprintf/
// http://en.wikipedia.org/wiki/Mantel_test
// http://www.bio.umontreal.ca/legendre/reprints/Partial_Mantel_paper.pdf
//
//


