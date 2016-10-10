#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/math_bot.h"
#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/solver_bot/sparse/sparse.h"

#include "sparse_struct.h"

#define POW2(a) ((a)*(a))

int nb_sparse_spy_plot_as_png(const nb_sparse_t *const A,
			       const char* url, uint32_t size,
			       bool enable_zeros_allocated,
			       bool enable_color)
{
//	/* Spyplot properties */
//	int border = 3;
//	int rgb_border[] = {120, 120, 120};
//	int s_palete = 4;
//	float p_palete[] = {0, 0.0001, 0.05, 0.2};
//	int rgb_palete[12];
//	if(enable_color){
//		int rgb [] =
//			{0, 0, 0,        /* Black */
//			 99, 0, 220,     /* Purple */
//			 255, 60, 0,     /* Orange */
//			 255, 255, 0};   /* Yellow */
//		memcpy(rgb_palete, rgb, 12 * sizeof(int));
//	}else{
//		int rgb[] = 
//			{255, 255, 255,    /* White */
//			 200, 200, 200,    /* Dark gray */
//			 50, 50, 50,       /* Gray */
//			 0, 0, 0};         /* Black */
//		memcpy(rgb_palete, rgb, 12 * sizeof(int));
//	}
//
//	/* Create spyplot */
//	int i, j, u, v;                   /* Iterative variables */
//	int plot_size = size-2*border;
//	float cell_size = ((float)A->N)/plot_size;
//	float cell_max = 0;
//	float* cells = nb_allocate_zero_mem(POW2(plot_size) * sizeof(float));
//	for(i=0; i<plot_size; i++){
//		for(j=0; j<plot_size; j++){
//			/* Start iterations over matrix */
//			for(u=floor(i*cell_size); u<ceil((i+1)*cell_size) && u < A->N; u++){
//				v = 0;
//				while(v < A->rows_size[u] && 
//				      A->rows_index[u][v] < ceil((j+1)*cell_size)){
//					if(!enable_zeros_allocated){
//						if(A->rows_values[u][v] == 0.0){
//							v++;
//							continue;
//						}
//					}
//					if(A->rows_index[u][v] >= floor(j*cell_size)){
//						float value = 1;
//						if(A->N > plot_size){
//							if(u == floor(i*cell_size))
//								value *= ceil(i*cell_size)-i*cell_size;
//							else if(u == ceil((i+1)*cell_size))
//								value *= (i+1)*cell_size-floor((i+1)*cell_size);
//	    
//							if(A->rows_index[u][v] == floor(j*cell_size))
//								value *= ceil(j*cell_size)-j*cell_size;
//							else if(A->rows_index[u][v] == ceil((j+1)*cell_size)-1)
//								value *= (j+1)*cell_size-floor((j+1)*cell_size);
//						}
//						cells[i*plot_size+j] += value;	      
//					}
//					v++;
//				}
//			}
//			/* Finish iterations over matrix */
//			if(cells[i*plot_size+j] > cell_max)
//				cell_max = cells[i*plot_size+j];
//		}      
//	}
//
//	/* Normalize cells matrix */
//	for(i=0; i<plot_size*plot_size; i++)
//		cells[i] /= cell_max;
//
//	/* Declare structures and variables to be used */
//	png_structp png_ptr = NULL;
//	png_infop info_ptr = NULL;
//	png_byte** row_pointers = NULL;
//
//	/* Image Properties */
//	int pixel_size = 3;
//	int depth = 8;
//	uint w = size;
//	uint h = size;
//  
//	/* Open url to write PNG file */
//	FILE *fp = fopen(url, "wb");
//
//	/* Verify all pointers are working */
//	if(fp == NULL){
//		printf("ERROR: Open url to save PNG failed\n");
//		return 1;
//	}
//  
//	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//	if(png_ptr == NULL){
//		printf("ERROR: PNG create write struct failed.\n");
//		fclose(fp);
//		return 2;
//	}
//  
//	info_ptr = png_create_info_struct(png_ptr);
//	if(info_ptr == NULL){
//		printf("ERROR: PNG create info struct failed.\n");
//		png_destroy_write_struct(&png_ptr, &info_ptr);
//		fclose(fp);
//		return 3;
//	}
//  
//	if(setjmp(png_jmpbuf(png_ptr))){
//		printf("ERROR: PNG failure.\n");
//		png_destroy_write_struct(&png_ptr, &info_ptr);
//		fclose(fp);
//		return 4;
//	}
//
//	/* Set image attributes */
//	png_set_IHDR(png_ptr, info_ptr, w, h, depth, 
//		     PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
//		     PNG_COMPRESSION_TYPE_DEFAULT,
//		     PNG_FILTER_TYPE_DEFAULT);
//
//	/* Initialize pointers */
//	uint8_t r, g, b;
//	row_pointers = png_nb_allocate_mem(png_ptr, h*sizeof(png_byte*));
//	for(i=0; i<h; i++){
//		png_byte *row=
//			png_nb_allocate_mem(png_ptr, sizeof(uint8_t)*w*pixel_size);
//		row_pointers[i] = row;
//		for(j=0; j<w; j++){
//			/* Create bitmap */
//			if(j < border || j > size-border-1 ||
//			   i < border || i > size-border-1){
//				/* Create matrix borders */
//				r = rgb_border[0];
//				g = rgb_border[1];
//				b = rgb_border[2];
//			}else{
//				int x = j-border;
//				int y = i-border;
//				float val = cells[y*plot_size+x];
//				int p = 0;
//				while(p < s_palete && val > p_palete[p])
//					p++;
//				if(p == 0){
//					r = rgb_palete[0];
//					g = rgb_palete[1];
//					b = rgb_palete[2];
//				}else if(p == s_palete){
//					r = rgb_palete[(p-1)*3];
//					g = rgb_palete[(p-1)*3+1];
//					b = rgb_palete[(p-1)*3+2];
//				}else{
//					float val_scaled = (val-p_palete[p-1])/(p_palete[p]-p_palete[p-1]);
//					r = (rgb_palete[p*3]-rgb_palete[(p-1)*3])*val_scaled 
//						+ rgb_palete[(p-1)*3];
//					g = (rgb_palete[p*3+1]-rgb_palete[(p-1)*3+1])*val_scaled
//						+ rgb_palete[(p-1)*3+1];
//					b = (rgb_palete[p*3+2]-rgb_palete[(p-1)*3+2])*val_scaled
//						+ rgb_palete[(p-1)*3+2];
//				}
//			}
//			/* End bitmap creation */
//			*row++ = r;
//			*row++ = g;
//			*row++ = b;
//		}
//	}
//
//	/* Write the image data */
//	png_init_io(png_ptr, fp);
//	png_set_rows(png_ptr, info_ptr, row_pointers);
//	png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
//
//	/* Free memory */
//	for(i=0; i<h; i++)
//		png_nb_free_mem(png_ptr, row_pointers[i]);
//	png_nb_free_mem(png_ptr, row_pointers);
//	png_destroy_write_struct(&png_ptr, &info_ptr);
//
//	nb_free_mem(cells);
//
//	/* Close url */
//	fclose(fp);
	return 0;
}
