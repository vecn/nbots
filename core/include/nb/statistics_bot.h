/**
 * @file statistics_bot.h
 * @brief Simple statistics routines.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 * @date 10 August 2015
 */

#ifndef __NB_STATISTICS_BOT_H__
#define __NB_STATISTICS_BOT_H__

#include <stdint.h>
	
uint32_t vcn_statistics_get_seed(void);
uint32_t vcn_statistics_lcg(uint32_t seed);
void vcn_statistics_random_permutation(uint32_t N, void *base, 
				       uint16_t type_size);
void vcn_statistics_runif(int n, double min, double max,
			  double *const out, uint64_t *const seed);

void vcn_statistics_rnorm(int n, double mean, double var,
			  double *const out, uint64_t *const seed1,
			  uint64_t *const seed2);

#endif
