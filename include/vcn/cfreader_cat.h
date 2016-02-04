/**
 * @file cfreader_cat.h
 * @brief Util to read customized file formats in plain text.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date 14 August 2015
 */

#ifndef __VCN_CFREADER_H__
#define __VCN_CFREADER_H__

#include <stdint.h>
#include <stdbool.h>

typedef struct vcn_cfreader_s vcn_cfreader_t;
vcn_cfreader_t* vcn_cfreader_create(const char* filename,
				    const char* line_comment_token);
char vcn_cfreader_read_int(vcn_cfreader_t *cfreader, int *val);
char vcn_cfreader_read_uint32_t(vcn_cfreader_t *cfreader, uint32_t *val);
char vcn_cfreader_read_float(vcn_cfreader_t *cfreader, float *val);
char vcn_cfreader_read_double(vcn_cfreader_t *cfreader, double *val);
char vcn_cfreader_read_bool(vcn_cfreader_t *cfreader, bool *val);
char* vcn_cfreader_read_and_allocate_string(vcn_cfreader_t *cfreader);
void vcn_cfreader_destroy(vcn_cfreader_t *cfreader);

#endif
