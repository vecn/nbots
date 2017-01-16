#ifndef __NB_CFREADER_BOT_H__
#define __NB_CFREADER_BOT_H__

#include <stdbool.h>
#include <stdint.h>

typedef struct nb_cfreader_s nb_cfreader_t;

uint32_t nb_cfreader_get_memsize(void);
void nb_cfreader_init(nb_cfreader_t *cfr);
nb_cfreader_t* nb_cfreader_create(void);
int nb_cfreader_open_file(nb_cfreader_t *cfr, const char* filename);
int nb_cfreader_add_line_comment_token(nb_cfreader_t *cfr, const char *token);
int nb_cfreader_read_int(nb_cfreader_t *cfr, int *val);
int nb_cfreader_read_uint(nb_cfreader_t *cfr, unsigned int *val);
int nb_cfreader_read_float(nb_cfreader_t *cfr, float *val);
int nb_cfreader_read_double(nb_cfreader_t *cfr, double *val);
int nb_cfreader_read_bool(nb_cfreader_t *cfr, bool *val);
char* nb_cfreader_read_and_allocate_string(nb_cfreader_t *cfr);
int nb_cfreader_read_var(nb_cfreader_t *cfr, char *var_name);
bool nb_cfreader_check_line(nb_cfreader_t *cfr, const char *line);
bool nb_cfreader_check_token(nb_cfreader_t *cfr, const char *token);
void nb_cfreader_close_file(nb_cfreader_t *cfr);
void nb_cfreader_finish(nb_cfreader_t *cfr);
void nb_cfreader_destroy(nb_cfreader_t *cfr);

#endif
