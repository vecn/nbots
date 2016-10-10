/**
 * @file cfreader_cat.h
 * @brief Util to read customized file formats in plain text.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n <a href="https://twitter.com/victore_cardoso"> @@victore_cardoso </a>
 * @date 14 August 2015
 */

#ifndef __NB_CFREADER_BOT_H__
#define __NB_CFREADER_BOT_H__

#include <stdbool.h>

typedef struct nb_cfreader_s nb_cfreader_t;
nb_cfreader_t* nb_cfreader_create(const char* filename,
				    const char* line_comment_token);
char nb_cfreader_read_int(nb_cfreader_t *cfreader, int *val);
char nb_cfreader_read_uint(nb_cfreader_t *cfreader, unsigned int *val);
char nb_cfreader_read_float(nb_cfreader_t *cfreader, float *val);
char nb_cfreader_read_double(nb_cfreader_t *cfreader, double *val);
char nb_cfreader_read_bool(nb_cfreader_t *cfreader, bool *val);
char* nb_cfreader_read_and_allocate_string(nb_cfreader_t *cfreader);
void nb_cfreader_destroy(nb_cfreader_t *cfreader);

#endif
