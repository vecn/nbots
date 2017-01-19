#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "nb/memory_bot.h"

#include "nb/io_bot.h"

#define LINE_LENGTH 256
#define MAX_LINE_COMMENT_TOKENS 5
#define MAX_ASSIGNMENT_TOKENS 5
#define READ_VAR_TOKEN_MAX_SIZE 20

static void init_bool_tokens(nb_cfreader_t *cfr);
static void init_string_tokens(nb_cfreader_t *cfr);
static int get_next_line(nb_cfreader_t* cfr);
static void set_end_of_line(nb_cfreader_t* cfr);
static void set_end_of_line_if_comment(nb_cfreader_t* cfr);
static void trim_start_of_line(nb_cfreader_t* cfr);
static int cast_bool(nb_cfreader_t *cfr, const char *string, bool *val);
static int get_string_value(nb_cfreader_t *cfr, char *line,
			    char *string);
static int trim_end_of_string(char *string, int size);

struct nb_cfreader_s {
	/* OPPORTUNITY: Support several tokens for comments, multiline comments,
	 * open_close_string (of several chars,
	 * e.g. open/close string with << >>), etc.
	 */
	uint8_t N_lct;                 /* Number of line comment tokens */
	uint8_t N_assign;                /* Number of assginment tokens */
	const char *lc_token[MAX_LINE_COMMENT_TOKENS];
	const char *assign_token[MAX_ASSIGNMENT_TOKENS];
	const char *bool_tokens[8];
	const char *string_tokens[4];
	FILE *fp;
	char buffer[LINE_LENGTH];
	char* line;
	uint32_t line_counter;
};

uint32_t nb_cfreader_get_memsize(void)
{
	return sizeof(nb_cfreader_t);
}

void nb_cfreader_init(nb_cfreader_t *cfr)
{
	cfr->N_lct = 0;
	cfr->N_assign = 0;
	init_bool_tokens(cfr);
	init_string_tokens(cfr);
	cfr->fp = NULL;
	cfr->buffer[0] = '\0';
	cfr->line = cfr->buffer;
	cfr->line_counter = 0;
}

static void init_bool_tokens(nb_cfreader_t *cfr)
{
	cfr->bool_tokens[0] = "TRUE";
	cfr->bool_tokens[1] = "FALSE";
	cfr->bool_tokens[2] = "true";
	cfr->bool_tokens[3] = "false";
	cfr->bool_tokens[4] = "T";
	cfr->bool_tokens[5] = "F";
	cfr->bool_tokens[6] = "1";
	cfr->bool_tokens[7] = "0";
}

static void init_string_tokens(nb_cfreader_t *cfr)
{
	cfr->string_tokens[0] = "\"";
	cfr->string_tokens[1] = "\"";
	cfr->string_tokens[2] = "'";
	cfr->string_tokens[3] = "'";
}

nb_cfreader_t* nb_cfreader_create(void)
{
	uint32_t memsize = nb_cfreader_get_memsize();
	nb_cfreader_t *cfr = nb_allocate_mem(memsize);
	nb_cfreader_init(cfr);
	return cfr;
}

int nb_cfreader_open_file(nb_cfreader_t *cfr, const char* filename)
{
	int status;
	cfr->fp = fopen(filename, "r");
	if (cfr->fp == NULL) {
		status = NB_CFREADER_FILE_NOT_FOUND;
		goto EXIT;
	}
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_add_line_comment_token(nb_cfreader_t *cfr, const char *token)
{
	int status;
	if (cfr->N_lct >= MAX_LINE_COMMENT_TOKENS) {
		status = NB_CFREADER_MAX_TOKENS_REACHED;
		goto EXIT;
	}

	cfr->lc_token[cfr->N_lct] = token;
	cfr->N_lct += 1;
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_add_assignment_token(nb_cfreader_t *cfr, const char *token)
{
	int status;
	if (cfr->N_assign >= MAX_ASSIGNMENT_TOKENS) {
		status = NB_CFREADER_MAX_TOKENS_REACHED;
		goto EXIT;
	}

	cfr->assign_token[cfr->N_assign] = token;
	cfr->N_assign += 1;
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

static int get_next_line(nb_cfreader_t* cfr)
{
	int status;
	while (fgets(cfr->buffer, LINE_LENGTH, cfr->fp) != NULL) {
		cfr->line_counter += 1;
		set_end_of_line(cfr);
		set_end_of_line_if_comment(cfr);

		cfr->line = cfr->buffer;
		trim_start_of_line(cfr);

		if (strlen(cfr->line) > 0) {
			status = NB_CFREADER_SUCCESS;
			goto EXIT;
		}
	}
	status = NB_CFREADER_EOF;
EXIT:
	return status;
}

static void set_end_of_line(nb_cfreader_t* cfr)
{
	char* pch = strchr(cfr->buffer, '\n');
	if (pch != NULL)
		pch[0] = '\0';
}

static void set_end_of_line_if_comment(nb_cfreader_t* cfr)
{
	for (int i = 0; i < cfr->N_lct; i++) {
		char *pch = strstr(cfr->buffer, cfr->lc_token[i]);
		if (pch != NULL)
			pch[0] = '\0';
	}
}

static void trim_start_of_line(nb_cfreader_t* cfr)
{
	while (cfr->line[0] == ' ' || cfr->line[0] == '\t')
		cfr->line = &(cfr->line[1]);
}

int nb_cfreader_read_int(nb_cfreader_t* cfr, int* val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	int len;
	if (sscanf(cfr->line, "%i%n", val, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_uint(nb_cfreader_t *cfr, uint32_t *val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	int len;
	if (sscanf(cfr->line, "%u%n", val, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_float(nb_cfreader_t* cfr, float* val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	int len;
	if (sscanf(cfr->line, "%f%n", val, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}
int nb_cfreader_read_double(nb_cfreader_t* cfr, double* val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	int len;
	if (sscanf(cfr->line, "%lf%n", val, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_bool(nb_cfreader_t* cfr, bool* val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	char token[256];
	int len;
	if (sscanf(cfr->line, "%s%n", token, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	status = cast_bool(cfr, token, val);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

static int cast_bool(nb_cfreader_t *cfr, const char *string, bool *val)
{
	int status = NB_CFREADER_BAD_INPUT;
	for (int i = 0; i < 4; i++) {
		if (0 == strcmp(string, cfr->bool_tokens[i * 2])) {
			*val = true;
			status = NB_CFREADER_SUCCESS;
			break;
		}
		if (0 == strcmp(string, cfr->bool_tokens[i*2+1])) {
			*val = false;
			status = NB_CFREADER_SUCCESS;
			break;
		}
	}
	return status;
}

int nb_cfreader_read_string(nb_cfreader_t *cfr, char *val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (0 == strlen(cfr->line)) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	status = get_string_value(cfr, cfr->line, val);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;

}

static int get_string_value(nb_cfreader_t *cfr, char *line,
			    char *string)
{
	char *pch_open = NULL;
	int token_id = 0;
	for (int i = 0; i < 2; i++) {
		char *aux;
		aux = strstr(line, cfr->string_tokens[i * 2]);
		if (aux != NULL) {
			if (pch_open != NULL) {
				if (aux < pch_open) {
					token_id = i;
					pch_open = aux;
				}
			} else {
				token_id = i;
				pch_open = aux;
			}
		}
	}
	int status;
	if (NULL == pch_open) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	int len = strlen(cfr->string_tokens[token_id*2]);
	pch_open += len;

	char *pch_close = strstr(pch_open, cfr->string_tokens[token_id*2+1]);
	if (NULL == pch_close) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	
	pch_close[0] = '\0';

	strcpy(string, pch_open);
	
	cfr->line = pch_close + 1;
	trim_start_of_line(cfr);

	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_token(nb_cfreader_t *cfr, char *val)
{
	int status = 1;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	int len;
	if (sscanf(cfr->line, "%s%n", val, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_tuple(nb_cfreader_t *cfr, char *var, char *val)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (0 == cfr->N_assign) {
		status = NB_CFREADER_ASSIGNMENT_TOKEN_UNDEF;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	/* Identify assignment */
	char *pch = NULL;
	int assignment_length;
	for (int i = 0; i < cfr->N_assign; i++) {
		pch = strstr(cfr->line, cfr->assign_token[i]);
		if (NULL != pch) {
			assignment_length = strlen(cfr->assign_token[i]);
			break;
		}
	}
	if (NULL == pch) {
		status = NB_CFREADER_NO_ASSIGNMENT;
		goto EXIT;
	}
	

	/* Read variable name */
	int var_size = pch - cfr->line;
	memcpy(var, cfr->line, var_size);

	var_size = trim_end_of_string(var, var_size);

	if (0 == var_size) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	/* Jump assignment token */
	cfr->line = pch + assignment_length;	
	trim_start_of_line(cfr);

	/* Read value */
	int len;
	if (sscanf(cfr->line, "%s%n", val, &len) != 1) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	if (0 == strlen(val)) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}
	cfr->line += len;
	trim_start_of_line(cfr);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

static int trim_end_of_string(char *string, int size)
{
	while ((string[size - 1] == ' ' || string[size - 1] == '\t') &&
	       size > 0) {
		string[size - 1] = '\0';
		size -= 1;
	}
	return size;
}

int nb_cfreader_read_var_int(nb_cfreader_t *cfr, const char *var, int *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}

	char *pch;
	long int N = strtol(val_readed, &pch, 10);
	if (pch == val_readed) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	*val = N;
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_var_uint(nb_cfreader_t *cfr, const char *var,
			      uint32_t *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}

	char *pch;
	unsigned long int N = strtoul(val_readed, &pch, 10);
	if (pch == val_readed) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	*val = N;
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_var_float(nb_cfreader_t *cfr, const char *var, float *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}

	char *pch;
	float N = strtof(val_readed, &pch);
	if (pch == val_readed) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	*val = N;
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_var_double(nb_cfreader_t *cfr, const char *var, double *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (0 != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}

	char *pch;
	double N = strtod(val_readed, &pch);
	if (pch == val_readed) {
		status = NB_CFREADER_BAD_INPUT;
		goto EXIT;
	}

	*val = N;
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_var_bool(nb_cfreader_t *cfr, const char *var, bool *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}

	status = cast_bool(cfr, val_readed, val);
EXIT:
	return status;
}


int nb_cfreader_read_var_token(nb_cfreader_t *cfr, const char *var, char *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}
	
	strcpy(val, val_readed);
	status = NB_CFREADER_SUCCESS;
EXIT:
	return status;
}

int nb_cfreader_read_var_string(nb_cfreader_t *cfr, const char *var, char *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (NB_CFREADER_SUCCESS != status)
		goto EXIT;

	if (0 != strcmp(var_readed, var)) {
		status = NB_CFREADER_DISTINCT_VARNAME;
		goto EXIT;
	}
	status = get_string_value(cfr, val_readed, val);
EXIT:
	return status;
}

int nb_cfreader_check_line(nb_cfreader_t *cfr, const char *line)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	trim_end_of_string(cfr->line, strlen(cfr->line));

	int cmp = strcmp(cfr->line, line);

	if (0 == cmp) {
		cfr->line[0] = '\0';
		status = NB_CFREADER_SUCCESS;
	} else {
		status = NB_CFREADER_BAD_INPUT;
	}
EXIT:
	return status;
}

int nb_cfreader_check_token(nb_cfreader_t *cfr, const char *token)
{
	int status;
	if (cfr->fp == NULL) {
		status = NB_CFREADER_NO_FILE_OPENED;
		goto EXIT;
	}

	if (strlen(cfr->line) == 0) {
		status = get_next_line(cfr);
		if (NB_CFREADER_EOF == status)
			goto EXIT;
	}

	int cmp = strncmp(cfr->line, token, strlen(token));

	if (0 == cmp) {
		int len = strlen(token);
		cfr->line += len;
		trim_start_of_line(cfr);
		status = NB_CFREADER_SUCCESS;
	} else {
		status = NB_CFREADER_BAD_INPUT;
	}
EXIT:
	return status;
}

void nb_cfreader_close_file(nb_cfreader_t *cfr)
{
	if (cfr->fp != NULL) {
		fclose(cfr->fp);
		cfr->fp = NULL;
	}
}

void nb_cfreader_finish(nb_cfreader_t *cfr)
{
	nb_cfreader_close_file(cfr);
}

void nb_cfreader_destroy(nb_cfreader_t* cfr)
{
	nb_free_mem(cfr);
}

void nb_cfreader_get_error_message(int error, char **msg)
{
	switch (error) {
	case NB_CFREADER_SUCCESS:
		*msg = "Successful read";
		break;
	case NB_CFREADER_EOF:
		*msg = "End of file";
		break;
	case NB_CFREADER_BAD_INPUT:
		*msg = "Incorrect or unexpected input";
		break;
	case NB_CFREADER_NO_FILE_OPENED:
		*msg = "There is not an open file";
		break;
	case NB_CFREADER_FILE_NOT_FOUND:
		*msg = "File not found";
		break;
	case NB_CFREADER_MAX_TOKENS_REACHED:
		*msg = "Max tokens reached";
		break;
	case NB_CFREADER_ASSIGNMENT_TOKEN_UNDEF:
		*msg = "There are not assignment tokens defined";
		break;
	case NB_CFREADER_NO_ASSIGNMENT:
		*msg = "There is not an assignment";
		break;
	case NB_CFREADER_DISTINCT_VARNAME:
		*msg = "Distinct variable name";
		break;
	}

}
