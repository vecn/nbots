/******************************************************************************
 *   Visual Cat: Visualization utilities.                                     *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "vcn/visual_cat.h"

struct vcn_palette_s{
  /* The palette defines a serie of RGB colors to
   * colorize values in [0,1]
   *
   *  c1    c2       c3         c4  <- RGB colors
   *   |_____|________|__________|
   *   0    0.25     0.57        1  <- Tics
   */
  uint8_t ntics;   /* Number of tics */
  float* tics;   /* Sorted tics in [0,1] */
  uint8_t* rgb;    /* RGB Colors definition */
};

static vcn_palette_t* palette_get_rainbow();
static vcn_palette_t* palette_get_sunset();
static vcn_palette_t* palette_get_french();

/************* Public functions **************************/
vcn_palette_t* vcn_palette_create(){
  return (vcn_palette_t*)calloc(1, sizeof(vcn_palette_t));
}

vcn_palette_t* vcn_palette_create_preset(int palette_id)
{
  if (VCN_PALETTE_RAINBOW == palette_id)
    return palette_get_rainbow();
  
  if (VCN_PALETTE_SUNSET == palette_id)
    return palette_get_sunset();

  if (VCN_PALETTE_FRENCH == palette_id)
    return palette_get_french();

  return NULL;
}

void vcn_palette_destroy(vcn_palette_t* palette){
  if(palette->ntics > 0){
    free(palette->tics);
    free(palette->rgb);
  }
  free(palette);
}

void vcn_palette_clear(vcn_palette_t *palette){
  if(palette->ntics > 0){
    free(palette->tics);
    free(palette->rgb);
    palette->tics = NULL;
    palette->rgb = NULL;
  }
  palette->ntics = 0;
}

void vcn_palette_add_colour(vcn_palette_t* palette, float tic,
                      uint8_t r, uint8_t g, uint8_t b){
  if(tic < 0) tic = 0;
  if(tic > 1) tic = 1;
  if(palette->ntics == 0){
    /* Insert first color */
    palette->tics = (float*)malloc(sizeof(float));
    palette->tics[0] = tic;
    palette->rgb = (uint8_t*)malloc(3);
    palette->rgb[0] = r;
    palette->rgb[1] = g;
    palette->rgb[2] = b;
    palette->ntics = 1;
  }else{
    /* Create a new space */
    float* tics = (float*)malloc(palette->ntics*sizeof(float));
    uint8_t* rgb = (uint8_t*)malloc(palette->ntics*3);
    memcpy(tics, palette->tics, palette->ntics*sizeof(float));
    memcpy(rgb, palette->rgb, palette->ntics*3);
    free(palette->tics);
    free(palette->rgb);
    palette->ntics += 1;
    palette->tics = (float*)malloc(palette->ntics*sizeof(float));
    palette->rgb = (uint8_t*)malloc(palette->ntics*3);
    memcpy(palette->tics, tics, (palette->ntics-1)*sizeof(float));
    memcpy(palette->rgb, rgb, (palette->ntics-1)*3);
    free(rgb);
    free(tics);
    /* Insert new color */
    palette->tics[palette->ntics-1] = 2;
    for(uint32_t i=0; i<palette->ntics; i++){
      if(tic < palette->tics[i]){
        float aux1 = tic;
        tic = palette->tics[i];
        palette->tics[i] = aux1;
        uint8_t aux2[3] = {r, g, b};
        r = palette->rgb[i * 3];
        g = palette->rgb[i*3+1];
        b = palette->rgb[i*3+2];
        palette->rgb[i * 3] = aux2[0];
        palette->rgb[i*3+1] = aux2[1];
        palette->rgb[i*3+2] = aux2[2];
      }
    }
  }
}

void vcn_palette_get_colour(const vcn_palette_t *const palette,
			    float factor,
			    uint8_t rgb[3])
{
  if (factor <= palette->tics[0]) {
    memcpy(rgb, palette->rgb, 3);
  } else if (factor >= palette->tics[palette->ntics-1]) {
    memcpy(rgb, &(palette->rgb[(palette->ntics-1)*3]), 3);
  } else {
    uint32_t i = 1;
    while(factor > palette->tics[i])
      i++;
    float w1 = (palette->tics[i]-factor)/(palette->tics[i]-palette->tics[i-1]);
    float w2 = (factor-palette->tics[i-1])/(palette->tics[i]-palette->tics[i-1]);
    rgb[0] = w1 * palette->rgb[(i-1) * 3] + w2 * palette->rgb[i * 3];
    rgb[1] = w1 * palette->rgb[(i-1)*3+1] + w2 * palette->rgb[i*3+1];
    rgb[2] = w1 * palette->rgb[(i-1)*3+2] + w2 * palette->rgb[i*3+2];
  }
}

/***************** Private functions *********************************/

static vcn_palette_t* palette_get_rainbow()
{
  vcn_palette_t* palette = vcn_palette_create();
  vcn_palette_add_colour(palette, 0.00f,   0,   0, 128);
  vcn_palette_add_colour(palette, 0.10f,   0,   0, 255);
  vcn_palette_add_colour(palette, 0.20f,   0, 128, 255);
  vcn_palette_add_colour(palette, 0.37f,   0, 255, 255);
  vcn_palette_add_colour(palette, 0.50f,   0, 255,   0);
  vcn_palette_add_colour(palette, 0.63f, 255, 255,   0);
  vcn_palette_add_colour(palette, 0.80f, 255, 128,   0);
  vcn_palette_add_colour(palette, 0.90f, 255,   0,   0);
  vcn_palette_add_colour(palette, 1.00f, 100,   0,   0);
  return palette;
}

static vcn_palette_t* palette_get_sunset()
{
  vcn_palette_t* palette = vcn_palette_create();
  vcn_palette_add_colour(palette, 0.00f,   0,   0,   0);
  vcn_palette_add_colour(palette, 0.15f,  20,   0, 100);
  vcn_palette_add_colour(palette, 0.30f, 100,   0, 200);
  vcn_palette_add_colour(palette, 0.80f, 220, 100,   0);
  vcn_palette_add_colour(palette, 1.00f, 255, 255,   0);
  return palette;
}

static vcn_palette_t* palette_get_french()
{
  vcn_palette_t* palette = vcn_palette_create();
  vcn_palette_add_colour(palette, 0.00f,   0,   0, 150);
  vcn_palette_add_colour(palette, 0.20f,   0,   0, 255);
  vcn_palette_add_colour(palette, 0.30f, 180, 180, 255);
  vcn_palette_add_colour(palette, 0.50f, 255, 255, 255);
  vcn_palette_add_colour(palette, 0.70f, 255, 180, 180);
  vcn_palette_add_colour(palette, 0.80f, 255,   0,   0);
  vcn_palette_add_colour(palette, 1.00f, 150,   0,   0);
  return palette;
}
