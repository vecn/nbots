#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "nb/memory_bot.h"

#include "nb/cfreader_bot.h"

#define LINE_LENGTH 256
#define MAX_LINE_COMMENT_TOKENS 5
#define MAX_ASSIGNMENT_TOKENS 5
#define READ_VAR_TOKEN_MAX_SIZE 20

static int get_next_line(nb_cfreader_t* cfr);
static void set_end_of_line(nb_cfreader_t* cfr);
static void set_end_of_line_if_comment(nb_cfreader_t* cfr);
static void trim_start_of_line(nb_cfreader_t* cfr);
static void forward_pointer_in_the_line(nb_cfreader_t *cfr);
static int cast_bool(nb_cfreader_t *cfr, const char *string, bool *val);
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
	const char *bool_tokens[8];/* = {"TRUE", "FALSE", "true", "false",
				      "T", "F", "1", "0"};*/
	char open_string_token;
	char close_string_token;
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
	cfr->open_string_token  = '\"';
	cfr->close_string_token = '\"';
	cfr->fp = NULL;
	cfr->buffer[0] = '\0';
	cfr->line = cfr->buffer;
	cfr->line_counter = 0;
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
	int status = 1;
	cfr->fp = fopen(filename, "r");
	if (cfr->fp == NULL)
		goto EXIT;
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_add_line_comment_token(nb_cfreader_t *cfr, const char *token)
{
	int status = 1;
	if (cfr->N_lct >= MAX_LINE_COMMENT_TOKENS)
		goto EXIT;

	cfr->lc_token[cfr->N_lct] = token;
	cfr->N_lct += 1;
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_add_assignment_token(nb_cfreader_t *cfr, const char *token)
{
	int status = 1;
	if (cfr->N_assign >= MAX_ASSIGNMENT_TOKENS)
		goto EXIT;

	cfr->assign_token[cfr->N_assign] = token;
	cfr->N_assign += 1;
	status = 0;
EXIT:
	return status;
}

static int get_next_line(nb_cfreader_t* cfr)
{
	int status = 1;
	while (fgets(cfr->buffer, LINE_LENGTH, cfr->fp) != NULL) {
		cfr->line_counter += 1;
		set_end_of_line(cfr);
		set_end_of_line_if_comment(cfr);
		trim_start_of_line(cfr);

		if (strlen(cfr->line) > 0) {
			status = 0;
			break;
		}
	}
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
	cfr->line = cfr->buffer;
	while (cfr->line[0] == ' ' || cfr->line[0] == '\t')
		cfr->line = &(cfr->line[1]);
}

int nb_cfreader_read_int(nb_cfreader_t* cfr, int* val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	if (sscanf(cfr->line, "%i", val) != 1)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
	status = 0;
EXIT:
	return status;
}

static void forward_pointer_in_the_line(nb_cfreader_t *cfr)
{
	while (cfr->line[0] != ' ' && cfr->line[0] != '\t' &&
	       cfr->line[0] != '\0')
		cfr->line = &(cfr->line[1]);
	while (cfr->line[0] == ' ' || cfr->line[0] == '\t')
		cfr->line = &(cfr->line[1]);
}

int nb_cfreader_read_uint(nb_cfreader_t *cfr, uint32_t *val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	if (sscanf(cfr->line, "%u", val) != 1)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_float(nb_cfreader_t* cfr, float* val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	if (sscanf(cfr->line, "%f", val) != 1)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
	status = 0;
EXIT:
	return status;
}
int nb_cfreader_read_double(nb_cfreader_t* cfr, double* val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	if (sscanf(cfr->line, "%lf", val) != 1)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_bool(nb_cfreader_t* cfr, bool* val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if(strlen(cfr->line) == 0)
		get_next_line(cfr);

	char token[256];
	if (sscanf(cfr->line, "%s", token) != 1)
		goto EXIT;

	status = cast_bool(cfr, token, val);
	if (0 != status)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
EXIT:
	return status;
}

static int cast_bool(nb_cfreader_t *cfr, const char *string, bool *val)
{
	int status = 1;
	for (int i = 0; i < 4; i++) {
		if (0 == strcmp(string, cfr->bool_tokens[i * 2])) {
			*val = true;
			status = 0;
			break;
		}
		if (0 == strcmp(string, cfr->bool_tokens[i*2+1])) {
			*val = false;
			status = 0;
			break;
		}
	}
	return status;
}

char* nb_cfreader_read_and_allocate_string(nb_cfreader_t* cfr)
{
	char* string = NULL;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);
	/* Allocate and read string */
	bool throw_open_warning = false;
	if (cfr->line[0] != cfr->open_string_token)
		throw_open_warning = true;
	else
		cfr->line = &(cfr->line[1]);

	char* ptr = cfr->line;
	uint32_t char_counter = 0;
	bool throw_close_warning = false;
	while (ptr[0] != cfr->close_string_token) {
		if (ptr[0] == '\0') {
			throw_close_warning = true;
			break;
		}
		char_counter ++;
		ptr = &(ptr[1]);
	}

	if (char_counter > 0) {
		string = nb_allocate_mem(char_counter+1);
		memcpy(string, cfr->line, char_counter);
		string[char_counter] = '\0';
		cfr->line = &(cfr->line[char_counter]);
		/* Throw warnings */
		if (throw_open_warning && throw_close_warning) {
			printf("WARNING nb_cfreader: ");
			printf("Missing open/close tokens '%c/%c' of string '%s' at line %i.\n",
			       cfr->open_string_token,
			       cfr->close_string_token,
			       string,
			       cfr->line_counter);
		} else if (throw_open_warning) {
			printf("WARNING nb_cfreader: ");
			printf("Missing open token '%c' of string '%s' at line %i.\n",
			       cfr->open_string_token,
			       string,
			       cfr->line_counter);
		} else if (throw_close_warning) {
			printf("WARNING nb_cfreader: ");
			printf("Missing close token '%c' of string '%s' at line %i.\n",
			       cfr->close_string_token,
			       string,
			       cfr->line_counter);
		}
	}
	if (ptr[0] == cfr->close_string_token)
		cfr->line = &(cfr->line[1]);
  
	forward_pointer_in_the_line(cfr);

EXIT:
	return string;
}

int nb_cfreader_read_token(nb_cfreader_t *cfr, char *val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	if (sscanf(cfr->line, "%s", val) != 1)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_tuple(nb_cfreader_t *cfr, char *var, char *val)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (0 == cfr->N_assign)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

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
	if (NULL == pch)
		goto EXIT;
	

	/* Read variable name */
	int var_size = pch - cfr->line;
	memcpy(var, cfr->line, var_size);

	var_size = trim_end_of_string(var, var_size);

	if (0 == var_size)
		goto EXIT;

	/* Jump assignment token */
	cfr->line = pch + assignment_length;	
	forward_pointer_in_the_line(cfr);

	/* Read value */
	if (sscanf(cfr->line, "%s", val) != 1)
		goto EXIT;

	if (0 == strlen(val))
		goto EXIT;
	
	forward_pointer_in_the_line(cfr);

	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_var_int(nb_cfreader_t *cfr, const char *var, int *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (0 != status)
		goto EXIT;

	status = 1;
	if (0 != strcmp(var_readed, var))
		goto EXIT;

	char *pch;
	long int N = strtol(val_readed, &pch, 10);
	if (pch == val_readed)
		goto EXIT;

	*val = N;
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_var_uint(nb_cfreader_t *cfr, const char *var,
			      uint32_t *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (0 != status)
		goto EXIT;

	status = 1;
	if (0 != strcmp(var_readed, var))
		goto EXIT;

	char *pch;
	unsigned long int N = strtoul(val_readed, &pch, 10);
	if (pch == val_readed)
		goto EXIT;

	*val = N;
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_var_float(nb_cfreader_t *cfr, const char *var, float *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (0 != status)
		goto EXIT;

	status = 1;
	if (0 != strcmp(var_readed, var))
		goto EXIT;

	char *pch;
	float N = strtof(val_readed, &pch);
	if (pch == val_readed)
		goto EXIT;

	*val = N;
	status = 0;
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

	status = 1;
	if (0 != strcmp(var_readed, var))
		goto EXIT;

	char *pch;
	double N = strtod(val_readed, &pch);
	if (pch == val_readed)
		goto EXIT;

	*val = N;
	status = 0;
EXIT:
	return status;
}

int nb_cfreader_read_var_bool(nb_cfreader_t *cfr, const char *var, bool *val)
{
	char var_readed[READ_VAR_TOKEN_MAX_SIZE];
	char val_readed[READ_VAR_TOKEN_MAX_SIZE];
	int status = nb_cfreader_read_tuple(cfr, var_readed, val_readed);
	if (0 != status)
		goto EXIT;

	status = 1;
	if (0 != strcmp(var_readed, var))
		goto EXIT;

	status = cast_bool(cfr, val_readed, val);
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

bool nb_cfreader_check_line(nb_cfreader_t *cfr, const char *line)
{
	bool out = false;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	out = (0 == strcmp(cfr->line, line));

	if (out)
		get_next_line(cfr);
EXIT:
	return out;
}

bool nb_cfreader_check_token(nb_cfreader_t *cfr, const char *token)
{
	bool out = false;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	out = (0 == strcmp(cfr->line, token));

	if (out) {
		int len = strlen(token);
		cfr->line += len;
		forward_pointer_in_the_line(cfr);
	}
EXIT:
	return out;
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
