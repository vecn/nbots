/******************************************************************************
 *   Graph's Bot: Fast graph theory tools.                                    *
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
 * @file graph_bot-coming_soon.h
 * @brief Next features to be added to the Graph's Bot.
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
 *
 * @mainpage Graph's Bot
 * A graph operations for numerical analysis.
 */

#ifndef __VCN_GRAPH_BOT_COOMING_SOON_H__
#define __VCN_GRAPH_BOT_COOMING_SOON_H__

#ifdef __cplusplus
extern "C"{
#endif

  void vcn_label_minimum_bandwidth
  /* Gibbs, Poole and Stockmeyer 1976 */
              (uint N_nodes,
	       const uint *const N_connections, 
	       uint** connectivity_matrix,
	       uint* perm, /* Output */
	       uint* iperm /* Output */);

  void vcn_label_spectral
  /* Labeling using the Friedler vector (Second eigenvector of Q) */
              (uint N_nodes,
	       const uint *const N_connections, 
	       uint** connectivity_matrix,
	       uint* perm, /* Output */
	       uint* iperm /* Output */);

  void vcn_nested_dissection(uint N_nodes,
			     uint* N_connections, 
			     uint** connectivity_matrix,
                             uint* perm, uint* iperm);
  

#ifdef __cplusplus
}
#endif
