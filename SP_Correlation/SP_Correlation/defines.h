/**
 * @file defines.h
 * @author Bryant Adams
 * @date 3/20/12
 */

#ifndef SC_defines_h
#define SC_defines_h

/**
 * @brief To supply a 'random seed' that'll yield an identity permutation
 */
#define SEED_IDENTITY 0

/**
 * @brief To use the existing seed
 */
#define SEED_RANDOM -1

/**
 * @brief Decimal places to display for field/list printing
 */
#define DISPLAY_DIGITS 3

/**
 * @brief When converting Int arithmetic to Float arithmetic
 */
#define FLOATIFY 1.0

/**
 * @brief Speed up by only operating on upper diagonal
 */
#define UPPER_ONLY false
#endif
