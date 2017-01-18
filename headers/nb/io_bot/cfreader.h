#ifndef __NB_IO_BOT_CFREADER_H__
#define __NB_IO_BOT_CFREADER_H__

#include <stdbool.h>
#include <stdint.h>

typedef struct nb_cfreader_s nb_cfreader_t;

enum {
	NB_CFREADER_SUCCESS,
	NB_CFREADER_EOF,
	NB_CFREADER_BAD_INPUT,
	NB_CFREADER_NO_FILE_OPENED,
	NB_CFREADER_FILE_NOT_FOUND,
	NB_CFREADER_MAX_TOKENS_REACHED,
	NB_CFREADER_ASSIGNMENT_TOKEN_UNDEF,
	NB_CFREADER_NO_ASSIGNMENT,
	NB_CFREADER_DISTINCT_VARNAME
};

uint32_t nb_cfreader_get_memsize(void);
void nb_cfreader_init(nb_cfreader_t *cfr);
nb_cfreader_t* nb_cfreader_create(void);
int nb_cfreader_open_file(nb_cfreader_t *cfr, const char* filename);
int nb_cfreader_add_line_comment_token(nb_cfreader_t *cfr, const char *token);
int nb_cfreader_add_assignment_token(nb_cfreader_t *cfr, const char *token);
int nb_cfreader_read_int(nb_cfreader_t *cfr, int *val);
int nb_cfreader_read_uint(nb_cfreader_t *cfr, uint32_t *val);
int nb_cfreader_read_float(nb_cfreader_t *cfr, float *val);
int nb_cfreader_read_double(nb_cfreader_t *cfr, double *val);
int nb_cfreader_read_bool(nb_cfreader_t *cfr, bool *val);
int nb_cfreader_read_string(nb_cfreader_t *cfr, char *val);

int nb_cfreader_read_token(nb_cfreader_t *cfr, char *val);
int nb_cfreader_read_tuple(nb_cfreader_t *cfr, char *var, char *val);

int nb_cfreader_read_var_int(nb_cfreader_t *cfr, const char *var, int *val);
int nb_cfreader_read_var_uint(nb_cfreader_t *cfr, const char *var,
			      uint32_t *val);
int nb_cfreader_read_var_float(nb_cfreader_t *cfr, const char *var, float *val);
int nb_cfreader_read_var_double(nb_cfreader_t *cfr, const char *var,
				double *val);
int nb_cfreader_read_var_bool(nb_cfreader_t *cfr, const char *var, bool *val);
int nb_cfreader_read_var_token(nb_cfreader_t *cfr, const char *var, char *val);
int nb_cfreader_read_var_string(nb_cfreader_t *cfr, const char *var,
				char *val);

int nb_cfreader_check_line(nb_cfreader_t *cfr, const char *line);
int nb_cfreader_check_token(nb_cfreader_t *cfr, const char *token);
void nb_cfreader_close_file(nb_cfreader_t *cfr);
void nb_cfreader_finish(nb_cfreader_t *cfr);
void nb_cfreader_destroy(nb_cfreader_t *cfr);

void nb_cfreader_get_error_message(int error, char **msg);

#endif
