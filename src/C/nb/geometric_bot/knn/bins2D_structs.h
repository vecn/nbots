#ifndef __NB_GEOMETRIC_BOT_KNN_BIN2D_STRUCT_H__
#define __NB_GEOMETRIC_BOT_KNN_BIN2D_STRUCT_H__

#include <stdint.h>
#include "nb/container_bot/container.h"

typedef struct {
	int32_t x, y;
	vcn_container_t* points;
} bin2D_t;

struct vcn_bins2D_s {
	/* Structure to sort vertices
	 *
	 *              GRID (Spatial structure)
	 *   ___________________________________________
	 *  |  . |  . |    | .  #### .: |    |   .|    |
	 *  |____|____|__._|____####____|___ |____|____|
	 *  |    |    |    |  . #### .  |  . | ...|    |
	 *  |___.|____|____|____####.___|____|____|____|
	 *  |.   |    |   .|    ####    |    |    |  <-+--- Grid Cells 
	 *  |__._|____|____|____####____|___ |____|____|
	 *  |   .|    |  . |   .####  . |    |    |.   |
	 *  |    |    |    |    #### 0  |    |    |    |    Shifted cells to
	 *  |##########################################| <- avoid double zero.
	 *  |##########################################|    
	 *  |  . |    |    | .  #### .  | .  |    | .<-+---- Points
	 *  |____|____|____|__._####____|___ |____|____|   /
	 *  |    |  . |    |    ####    |.   |    |    |  /
	 *  |____|____|____|____####____|____|____|____| /
	 *  |    |   .|    |   .####    |    |    |    |/
	 *  |____|____|____|____####____|___ |____|____/
	 *  |.   |    |    |    #### .  |    |    |  ./|
	 *  |____|____|____|__._####____|____|____|____|
	 *           
	 *       |<-->|
	 *        Size
	 */
	double size_of_bins;
	uint32_t length;
	vcn_container_t* bins;
	bool (*filter)(const vcn_point2D_t *const p_ref, 
		       const vcn_point2D_t *const p,
		       const void *const data);
	const void *filter_data;
	void (*destroy)(void*);
	void (*destroy_attr)(void*);
};

#endif
