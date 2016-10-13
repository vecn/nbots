#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/image_bot.h"
#include "nb/graphics_bot.h"
#include "nb/solver_bot/sparse/sparse.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))

int nb_sparse_spy_plot_as_png(const nb_sparse_t *const A,
			      const char* url, uint32_t size,
			      bool enable_zeros_allocated,
			      bool enable_color)
{
	/* Spyplot properties */
	int border = 3;
	
	int plot_size = size-2*border;
	float cell_size = ((float)A->N)/plot_size;
	float cell_max = 0;
	float* cells = nb_allocate_zero_mem(POW2(plot_size) * sizeof(float));
	create_spyplot();
	
	/* Normalize cells matrix */
	for(i=0; i<plot_size*plot_size; i++)
		cells[i] /= cell_max;
	
	
	/* Initialize pointers */
	uint8_t r, g, b;
	row_pointers = png_nb_allocate_mem(png_ptr, h*sizeof(png_byte*));
	for(i=0; i<h; i++){
		png_byte *row=
			png_nb_allocate_mem(png_ptr, sizeof(uint8_t)*w*pixel_size);
		row_pointers[i] = row;
		for(j=0; j<w; j++){
			/* Create bitmap */
			if(j < border || j > size-border-1 ||
			   i < border || i > size-border-1){
				/* Create matrix borders */
				r = rgb_border[0];
				g = rgb_border[1];
				b = rgb_border[2];
			}else{
				int x = j-border;
				int y = i-border;
				float val = cells[y*plot_size+x];
				int p = 0;
				while(p < s_palete && val > p_palete[p])
					p++;
				if(p == 0){
					r = rgb_palete[0];
					g = rgb_palete[1];
					b = rgb_palete[2];
				}else if(p == s_palete){
					r = rgb_palete[(p-1)*3];
					g = rgb_palete[(p-1)*3+1];
					b = rgb_palete[(p-1)*3+2];
				}else{
					float val_scaled = (val-p_palete[p-1])/(p_palete[p]-p_palete[p-1]);
					r = (rgb_palete[p*3]-rgb_palete[(p-1)*3])*val_scaled 
						+ rgb_palete[(p-1)*3];
					g = (rgb_palete[p*3+1]-rgb_palete[(p-1)*3+1])*val_scaled
						+ rgb_palete[(p-1)*3+1];
					b = (rgb_palete[p*3+2]-rgb_palete[(p-1)*3+2])*val_scaled
						+ rgb_palete[(p-1)*3+2];
				}
			}
			/* End bitmap creation */
			*row++ = r;
			*row++ = g;
			*row++ = b;
		}
	}
	
	/* Write the image data */
	png_init_io(png_ptr, fp);
	png_set_rows(png_ptr, info_ptr, row_pointers);
	png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
	
	/* Free memory */
	for(i=0; i<h; i++)
		png_nb_free_mem(png_ptr, row_pointers[i]);
	png_nb_free_mem(png_ptr, row_pointers);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	
	nb_free_mem(cells);
	
	/* Close url */
	fclose(fp);
	return 0;
}

static void create_spyplot()
{
	for(i=0; i<plot_size; i++){
		for(j=0; j<plot_size; j++){
			/* Start iterations over matrix */
			for(u=floor(i*cell_size); u<ceil((i+1)*cell_size) && u < A->N; u++){
				v = 0;
				while(v < A->rows_size[u] && 
				      A->rows_index[u][v] < ceil((j+1)*cell_size)){
					if(!enable_zeros_allocated){
						if(A->rows_values[u][v] == 0.0){
							v++;
							continue;
						}
					}
					if(A->rows_index[u][v] >= floor(j*cell_size)){
						float value = 1;
						if(A->N > plot_size){
							if(u == floor(i*cell_size))
								value *= ceil(i*cell_size)-i*cell_size;
							else if(u == ceil((i+1)*cell_size))
								value *= (i+1)*cell_size-floor((i+1)*cell_size);
		    
							if(A->rows_index[u][v] == floor(j*cell_size))
								value *= ceil(j*cell_size)-j*cell_size;
							else if(A->rows_index[u][v] == ceil((j+1)*cell_size)-1)
								value *= (j+1)*cell_size-floor((j+1)*cell_size);
						}
						cells[i*plot_size+j] += value;	      
					}
					v++;
				}
			}
			/* Finish iterations over matrix */
			if(cells[i*plot_size+j] > cell_max)
				cell_max = cells[i*plot_size+j];
		}      
	}
}
