/**
 * @file array.h
 * @brief Simple math routines.
 * @author Victor Eduardo Cardoso Nungaray
 * @n @@victore_cardoso
 */

#ifndef __NB_CONTAINER_BOT_ARRAY_H__
#define __NB_CONTAINER_BOT_ARRAY_H__

#include <stdint.h>

void* nb_array_get(const void *const base, uint16_t type_size,
		    uint32_t i);
void nb_swap(void* base, uint32_t i, uint32_t j, uint16_t type_size);
void nb_qsort(void *base, uint32_t N, uint16_t type_size,
	       int8_t (*compare)(const void *const,
				 const void *const));
void nb_qsort_wd(void *base, uint32_t N, uint16_t type_size,
		  int8_t (*compare)(const void *const,
				    const void *const,
				    const void *const data),
		  const void *const data);
uint32_t nb_array_get_min_id(const void *const array, uint32_t N,
			      uint16_t type_size,
			      int8_t (*compare)(const void *const,
						const void *const));
uint32_t nb_array_get_min_id_wd(const void *const array,
				 uint32_t N,
				 uint16_t type_size,
				 int8_t (*compare)(const void *const,
						   const void *const,
						   const void *const d),
				 const void *const data);
uint32_t nb_array_get_max_id(const void *const array, uint32_t N,
			      uint16_t type_size,
			      int8_t (*compare)(const void *const,
						const void *const));
uint32_t nb_array_get_max_id_wd(const void *const array, uint32_t N,
				 uint16_t type_size,
				 int8_t (*compare)(const void *const,
						   const void *const,
						   const void *const d),
				 const void *const data);
void nb_array_get_min_max_ids(const void *const array, uint32_t N,
			       uint16_t type_size,
			       int8_t (*compare)(const void *const,
						 const void *const),
			       uint32_t *min, uint32_t *max);
void nb_array_get_min_max_ids_wd(const void *const array, uint32_t N,
				  uint16_t type_size,
				  int8_t (*compare)(const void *const,
						    const void *const,
						    const void *const d),
				  const void *const data,
				  uint32_t *min, uint32_t *max);
int8_t nb_compare_char(const void *const a, const void *const b);
int8_t nb_compare_float(const void *const a, const void *const b);
int8_t nb_compare_double(const void *const a, const void *const b);
int8_t nb_compare_int8(const void *const a, const void *const b);
int8_t nb_compare_int16(const void *const a, const void *const b);
int8_t nb_compare_int32(const void *const a, const void *const b);
int8_t nb_compare_int64(const void *const a, const void *const b);
int8_t nb_compare_uint8(const void *const a, const void *const b);
int8_t nb_compare_uint16(const void *const a, const void *const b);
int8_t nb_compare_uint32(const void *const a, const void *const b);
int8_t nb_compare_uint64(const void *const a, const void *const b);

#endif
