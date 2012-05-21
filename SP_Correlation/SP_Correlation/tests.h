/**
 * @file tests.h
 * @author Bryant Adams
 * @date 3/20/12
 */

#ifndef SC_tests_h
#define SC_tests_h
#include <stdbool.h>
/**
 * @brief Create a number of data files to test correlation with
 * @sideeffect Creates Test1a.tdv, Test1b.tdv, Test2a.tdv, etc
 */
void makeTestFiles(void);

/**
 * @brief Test functions for correct behavior
 */
bool runTests(void);

void reportStart(char* function);
bool reportEnd(bool good, char* explanation);
bool reportFake(void);

/**
 * @brief Generate random number in a given range
 */
bool testRandInRange(void);

/**
 * @brief Check if inputs are in desired order
 */
bool testInOrder(void);

/**
 * @brief Exchange contents of two float variables
 */
bool testSwapF(void);

/**
 * @brief Exchange contents of two int variables
 */
bool testSwapI(void);

#pragma mark Correlation

/**
 * @brief Initialize a correlation aggregate
 */
bool testInitializeCA(void);

/**
 * @brief Add data to a correlation aggregate
 */
bool testAugmentCAByValues(void);

/**
 * @brief Extract result from correlation aggregate
 */
bool testFinishCorrelation(void);

#pragma mark Permutation

/**
 * @brief Creates a permutation of a given size
 */
bool testMakePerm(void);

/**
 * @brief Change a permutation
 */
bool testModifyPermPermutify(void);

#pragma mark Lists

/**
 * @brief Create a List by copying another list
 */
bool testMakeListFromList(void);

/**
 * @brief Display a list on screen
 */
bool  testDisplayList(void);

/**
 * @brief Output a list to file
 */
bool  testSaveListToTDV(void);

/**
 * @brief Center a list about its mean
 */
bool  testModifyListMeanify(void);

/**
 * @brief Sort a list (heapsort)
 */
bool  testModifyListSortify(void);

/**
 * @brief Utility function for sorting
 */
bool  testHeapSortSiftDown(void);


bool testComputeRankInList(void);


#pragma mark Fields

/**
 * @brief Create a symmetric Field with random data
 */
bool testMakeRandomField(void);

/**
 * @brief Load a tab-delimited-value file into a Field
 */
bool testMakeFieldFromTDV(void);

/**
 * @brief Save field data to file.
 */
bool  testSaveFieldToTDV(void);

/**
 * @brief Save field data to file.
 */
bool testSaveFieldToHTML(void);

/**
 * @brief Display field data on screen.
 */
bool testDisplayField(void);

/**
 * @brief Center a field about its mean
 */
bool testModifyFieldMeanify(void);

/**
 * @brief Convert field of data to field of ranks
 */
bool testModifyFieldRankify(void);

/**
 * @brief Augment a correlation aggregate with information from two fields
 */
bool testAugmentCAByFields(void);


#pragma mark Landscapes

bool testMakeLandscapeFromTDVs(void);

bool testModifyLandscapeMeanify(void);

bool testModifyLandscapeRankify(void);

#pragma mark Field->List
/**
 * @brief Flatten a 2D field to a 1D list
 */
bool testMakeListFromField(void);

bool testMakeListFromLandscape(void);

#pragma mark Correlation
/**
 * @brief Find correlation of two sets of fields
 */

bool testMantelR(void);
bool testMantelRPartial(void);

#pragma mark P value

bool testCorrelateAndFindP(void);

bool testCorrelatePartialAndFindP(void);

/**
 * @brief Creates files containing data on Spearman and Person correlation of inputs
 */
bool testProcessFilePairs(void);

bool testProcessFileTriples(void);

/**
 * @brief Saves statistical data to file
 */
bool testSaveData(void);


#endif
