/******************************************************************************
 *   FEM Bot Coming Soon: Next features to be added on Finit Elem Bot         *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 *                                                                            *
 *   License:                                                                 *
 *   This piece of code is in the PUBLIC DOMAIN, if your government does not  *
 *   recognizes such a dedication, then you are granted a perpetual and       *
 *   irrevocable license to copy and modify this file however you want. This  *
 *   does not imply any warranty.                                             *
 *   Attributions and feedback are always welcome.                            *
 ******************************************************************************/

/**
 * @file fem_bot-coming_soon.h
 * @brief Next features add to the Finite Element Bot.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date 10 August 2015
 *
 * @par License:@n
 * This piece of code is in the PUBLIC DOMAIN, if your government does not
 * recognizes such a dedication, then you are granted a perpetual and
 * irrevocable license to copy and modify this file however you want. This
 * does not imply any warranty.
 * @n Attributions and feedback are always welcome.
 */

#ifndef __VCN_FEM_BOT_COMING_SOON_H__
#define __VCN_FEM_BOT_COMING_SOON_H__

#ifdef __cplusplus
extern "C" {
#endif

  void vcn_fem_enrich_output
  (const char* vcn_dma_filename,
   const char* vcn_dma_filename_enhanced,
   const vcn_fem_elem_t *const elemtype,
   const vcn_fem_material_t *const material,
   bool add_main_strain,
   bool add_real_stress,
   bool add_main_stress,
   bool add_effective_stress,
   bool add_main_effective_stress,
   bool add_von_mises_stress,
   bool add_error_based_on_strain){
    /* Read DMA file */
    vcn_dma_t* dma = vcn_dma_create(vcn_dma_filename);
    /* Get data  (N_time X N_atr X N_elm) */
    
    /* Compute main strain */
    /* [PENDING...] */
    /* Compute real stress */
    /* [PENDING...] */
    /* Compute main stress */
    /* [PENDING...] */
    /* Compute effective stress */
    /* [PENDING...] */
    /* Compute main effective stress */
    /* [PENDING...] */
    /* Compute von mises stress */
    /* [PENDING...] */
    /* Compute error based on strain */
    /* [PENDING...] */
    /* Save DMA */
    
    /* Free DMA data */
    vcn_dma_destroy(dma);
  }
  
  void vcn_fem_analyze_reaction
  (const char* vcn_dma_filename,
   const char* data_filename,
   bool track_vertical_reaction,
   bool track_horizontal_reaction,
   bool track_reaction_magnitude,
   uint N_nodes, uint* nodes);
  
#ifdef __cplusplus
}
#endif

#endif
