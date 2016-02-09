/******************************************************************************
 *   Array: array utilities                                                   *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "vcn/container_bot/array.h"

#define MAX_CHUNK_SIZE 16 /* Bytes */

static int8_t swap_chunks(void* base, uint32_t i, uint32_t j,
			  uint32_t type_size, uint8_t swapped_bytes, 
			  uint8_t chunk_size);
static void swap_mem(void* base, uint32_t i, uint32_t j,
		     uint32_t type_size, uint8_t swapped_bytes,
		     uint8_t chunk_size);
static char* get_pos(const void *const base, uint32_t type_size,
		     uint32_t i);
static void swap_bytes(void* base, uint32_t i, uint32_t j,
		       uint32_t type_size, uint8_t swapped_bytes);
static void swap_byte(void* base, uint32_t i, uint32_t j,
		      uint32_t type_size, uint8_t swapped_bytes);
static void qsort_data(void *base, int32_t p, int32_t r, uint16_t type_size,
		       int8_t (*compare)(), const void *const data);
static int lomuto_partition(void *base, int32_t p, int32_t r, uint16_t type_size,
			    int8_t (*compare)(), const void *const data);
static uint32_t array_get_min(const void *const array, uint32_t N,
			      uint16_t type_size, int8_t (*compare)(),
			      const void *const data);
static uint32_t array_get_max(const void *const array, uint32_t N,
			      uint16_t type_size, int8_t (*compare)(),
			      const void *const data);
static void array_get_min_max(const void *const array, uint32_t N,
			      uint16_t type_size,
			      int8_t (*compare)(),
			      const void *const data,
			      uint32_t *min, uint32_t *max);

inline void* vcn_array_get(const void *const base, uint16_t type_size,
			   uint32_t i)
{
	return get_pos(base, type_size, i);
}

void vcn_swap(void* base, uint32_t i, uint32_t j, uint16_t type_size)
{
	int8_t swapped_bytes = swap_chunks(base, i, j, type_size, 0, 8);
	swapped_bytes = swap_chunks(base, i, j, type_size, swapped_bytes, 4);
	swapped_bytes = swap_chunks(base, i, j, type_size, swapped_bytes, 2);
	swap_bytes(base, i, j, type_size, swapped_bytes);
}

static inline int8_t swap_chunks(void* base, uint32_t i, uint32_t j,
				 uint32_t type_size, uint8_t swapped_bytes, 
				 uint8_t chunk_size)
{
	while (type_size - swapped_bytes > chunk_size) {
		swap_mem(base, i, j, type_size, swapped_bytes, chunk_size);
		swapped_bytes += chunk_size;
	}
	return swapped_bytes;
}

static inline void swap_mem(void* base, uint32_t i, uint32_t j,
			    uint32_t type_size, uint8_t swapped_bytes, 
			    uint8_t chunk_size)
{
	char aux[MAX_CHUNK_SIZE];
	register char *i_pos = get_pos(base, type_size, i) + swapped_bytes;
	register char *j_pos = get_pos(base, type_size, j) + swapped_bytes;
	memcpy(aux, i_pos, chunk_size);
	memcpy(i_pos, j_pos, chunk_size);
	memcpy(j_pos, aux, chunk_size);
}

static inline char* get_pos(const void *const base, uint32_t size, uint32_t i) 
{
	return ((char*)base + i * size);
}

static inline void swap_bytes(void* base, uint32_t i, uint32_t j,
			      uint32_t type_size, uint8_t swapped_bytes)
{
	while (swapped_bytes < type_size){
		swap_byte(base, i, j, type_size, swapped_bytes);
		swapped_bytes += 1;
	}
}

static inline void swap_byte(void* base, uint32_t i, uint32_t j,
			     uint32_t type_size, uint8_t swapped_bytes)
{
	register char *i_pos = get_pos(base, type_size, i) + swapped_bytes;
	register char *j_pos = get_pos(base, type_size, j) + swapped_bytes;
	register char aux = *(i_pos);
	*(i_pos) = *(j_pos);
	*(j_pos) = aux;

}

inline void vcn_qsort(void *base, uint32_t N, uint16_t type_size,
		      int8_t (*compare)(const void *const, const void *const))
{
	qsort_data(base, 0, N - 1, type_size, compare, NULL);
}

inline void vcn_qsort_wd(void *base, uint32_t N, uint16_t type_size,
			 int8_t (*compare)(const void *const,
					   const void *const,
					   const void *const data),
			 const void *const data)
{
	qsort_data(base, 0, N - 1, type_size, compare, data);
}

static void qsort_data(void *base, int32_t p, int32_t r, uint16_t type_size,
		       int8_t (*compare)(), const void *const data)
{
	if (p < r) {
		int k = lomuto_partition(base, p, r, type_size, compare, data);	
		qsort_data(base, p, k-1, type_size, compare, data);
		qsort_data(base, k+1, r, type_size, compare, data);
	}
}

static int lomuto_partition(void *base, int32_t p, int32_t r, uint16_t type_size,
			    int8_t (*compare)(), const void *const data)
{
	void *rptr = get_pos(base, type_size, r);
	int k = p;
	for (uint32_t j = p; j < r; j++) {
		void *jptr = get_pos(base, type_size, j);
		if (compare(jptr, rptr, data) < 0) {
			vcn_swap(base, k, j, type_size);
			k++;
		}
	}
	vcn_swap(base, k, r, type_size);
	return k;
}

inline uint32_t vcn_array_get_min_id(const void *const array, uint32_t N,
				     uint16_t type_size,
				     int8_t (*compare)(const void *const,
						       const void *const))
{
	return array_get_min(array, N, type_size, compare, NULL);
}

inline uint32_t vcn_array_get_min_id_wd(const void *const array,
					uint32_t N,
					uint16_t type_size,
					int8_t (*compare)(const void *const,
							  const void *const,
							  const void *const d),
					const void *const data)
{
	return array_get_min(array, N, type_size, compare, data);
	
}

static uint32_t array_get_min(const void *const array, uint32_t N,
			      uint16_t type_size, int8_t (*compare)(),
			      const void *const data)
{
	uint32_t min = 0;
	char *min_pos = (char*) array;
	for (uint32_t i = 1; i < N; i++) {
		char *pos = get_pos(array, type_size, i);
		if (compare(min_pos, pos, data) > 0) {
			min = i;
			min_pos = pos;
		}
	}
	return min;
}

inline uint32_t vcn_array_get_max_id(const void *const array, uint32_t N,
				     uint16_t type_size,
				     int8_t (*compare)(const void *const,
						       const void *const))
{
	return array_get_max(array, N, type_size, compare, NULL);
}

inline uint32_t vcn_array_get_max_id_wd(const void *const array,
					uint32_t N,
					uint16_t type_size,
					int8_t (*compare)(const void *const,
							  const void *const,
							  const void *const d),
					const void *const data)
{
	return array_get_max(array, N, type_size, compare, data);
}

static uint32_t array_get_max(const void *const array, uint32_t N,
			      uint16_t type_size, int8_t (*compare)(),
			      const void *const data)
{
	int32_t max = 0;
	char *max_pos = (char*) array;
	for (uint32_t i = 1; i < N; i++) {
		char *pos = get_pos(array, type_size, i);
		if (compare(pos, max_pos, data) > 0) {
			max = i;
			max_pos = pos;
		}
	}
	return max;

}

inline void vcn_array_get_min_max_ids(const void *const array, uint32_t N,
				      uint16_t type_size,
				      int8_t (*compare)(const void *const,
							const void *const),
				      uint32_t *min, uint32_t *max)
{
	array_get_min_max(array, N, type_size, compare, NULL, min, max);
}

inline void vcn_array_get_min_max_ids_wd(const void *const array, uint32_t N,
					 uint16_t type_size,
					 int8_t (*compare)(const void *const,
							   const void *const,
							   const void *const d),
					 const void *const data,
					 uint32_t *min, uint32_t *max)
{
	array_get_min_max(array, N, type_size, compare, data, min, max);
}

static void array_get_min_max(const void *const array, uint32_t N,
			      uint16_t type_size,
			      int8_t (*compare)(),
			      const void *const data,
			      uint32_t *min, uint32_t *max)
{
	*min = 0;
	*max = 0;
	char *min_pos = (char*) array;
	char *max_pos = (char*) array;
	for (uint32_t i = 1; i < N; i++) {
		char *pos = get_pos(array, type_size, i);
		if (compare(min_pos, pos, data) > 0) {
			*min = i;
			min_pos = pos;
		} else if (compare(pos, max_pos, data) > 0) {
			*max = i;
			max_pos = pos;
		}
	}
}

int8_t vcn_compare_char(const void *const restrict a, 
			const void *const restrict b)
{
	int8_t out = 0;
	if (*((char*)a) > *((char*)b))
		out = 1;
	else if (*((char*)a) < *((char*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_float(const void *const restrict a, 
			 const void *const restrict b)
{
	int8_t out = 0;
	if (*((float*)a) > *((float*)b))
		out = 1;
	else if (*((float*)a) < *((float*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_double(const void *const restrict a, 
			  const void *const restrict b)
{
	int8_t out = 0;
	if (*((double*)a) > *((double*)b))
		out = 1;
	else if (*((double*)a) < *((double*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_int8(const void *const restrict a, 
			const void *const restrict b)
{
	int8_t out = 0;
	if (*((int8_t*)a) > *((int8_t*)b))
		out = 1;
	else if (*((int8_t*)a) < *((int8_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_int16(const void *const restrict a, 
			 const void *const restrict b)
{
	int8_t out = 0;
	if (*((int16_t*)a) > *((int16_t*)b))
		out = 1;
	else if (*((int16_t*)a) < *((int16_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_int32(const void *const restrict a, 
			 const void *const restrict b)
{
	int32_t out = 0;
	if (*((int32_t*)a) > *((int32_t*)b))
		out = 1;
	else if (*((int32_t*)a) < *((int32_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_int64(const void *const restrict a, 
			 const void *const restrict b)
{
	int8_t out = 0;
	if (*((int64_t*)a) > *((int64_t*)b))
		out = 1;
	else if (*((int64_t*)a) < *((int64_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_uint8(const void *const restrict a, 
			 const void *const restrict b)
{
	int8_t out = 0;
	if (*((uint8_t*)a) > *((uint8_t*)b))
		out = 1;
	else if (*((uint8_t*)a) < *((uint8_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_uint16(const void *const restrict a, 
			  const void *const restrict b)
{
	int8_t out = 0;
	if (*((uint16_t*)a) > *((uint16_t*)b))
		out = 1;
	else if (*((uint16_t*)a) < *((uint16_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_uint32(const void *const restrict a, 
			  const void *const restrict b)
{
	int32_t out = 0;
	if (*((uint32_t*)a) > *((uint32_t*)b))
		out = 1;
	else if (*((uint32_t*)a) < *((uint32_t*)b))
		out = -1;
	return out;
}

int8_t vcn_compare_uint64(const void *const restrict a, 
			  const void *const restrict b)
{
	int8_t out = 0;
	if (*((uint64_t*)a) > *((uint64_t*)b))
		out = 1;
	else if (*((uint64_t*)a) < *((uint64_t*)b))
		out = -1;
	return out;
}
