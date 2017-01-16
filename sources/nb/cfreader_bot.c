#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "nb/memory_bot.h"

#include "nb/cfreader_bot.h"

#define LINE_LENGTH 256
#define MAX_LINE_COMMENT_TOKENS 5

static int get_next_line(nb_cfreader_t* cfr);
static void set_end_of_line(nb_cfreader_t* cfr);
static void set_end_of_line_if_comment(nb_cfreader_t* cfr);
static void trim_start_of_line(nb_cfreader_t* cfr);
static void forward_pointer_in_the_line(nb_cfreader_t *cfr);

struct nb_cfreader_s {
	/* OPPORTUNITY: Support several tokens for comments, multiline comments,
	 * open_close_string (of several chars,
	 * e.g. open/close string with << >>), etc.
	 */
	uint8_t N_lct;                     /* Number of line comment tokens */
     	                                             /* Line comment tokens */
	const char *lc_token[MAX_LINE_COMMENT_TOKENS];
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

int nb_cfreader_read_uint(nb_cfreader_t *cfr, unsigned int *val)
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

	char* ptr = cfr->line;
	uint32_t char_counter = 0;
	while (ptr[0] != ' ' && ptr[0] != '\t' && ptr[0] != '\0') {
		char_counter ++;
		ptr = &(ptr[1]);
	}
	char* value = nb_allocate_mem(char_counter+1);
	memcpy(value, cfr->line, char_counter);
	value[char_counter] = '\0';

	/* Read boolean value */
	bool is_a_correct_value = true;
	if (strcmp(value, "TRUE") == 0 || strcmp(value, "true") == 0 ||
	    strcmp(value, "T") == 0 || strcmp(value, "1") == 0)
		val[0] = true;
	else if (strcmp(value, "FALSE") == 0 ||
		 strcmp(value, "false") == 0 ||
		 strcmp(value, "F") == 0 ||
		 strcmp(value, "0") == 0)
		val[0] = false;
	else
		is_a_correct_value = false;

	forward_pointer_in_the_line(cfr);

	/* Show warning of wrong formating */
	if (!is_a_correct_value) {
		printf("WARNING Cfr: ");
		printf("The boolean at line %i has an invalid value: \"%s\").\n",
		       cfr->line_counter, value);
	}
	/* Free memory */
	nb_free_mem(value);

	status = 0;
EXIT:
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
			printf("WARNING Cfr: ");
			printf("Missing open/close tokens '%c/%c' of string '%s' at line %i.\n",
			       cfr->open_string_token,
			       cfr->close_string_token,
			       string,
			       cfr->line_counter);
		} else if (throw_open_warning) {
			printf("WARNING Cfr: ");
			printf("Missing open token '%c' of string '%s' at line %i.\n",
			       cfr->open_string_token,
			       string,
			       cfr->line_counter);
		} else if (throw_close_warning) {
			printf("WARNING Cfr: ");
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

int nb_cfreader_read_var(nb_cfreader_t *cfr, char *var_name)
{
	int status = 1;
	if (cfr->fp == NULL)
		goto EXIT;

	if (strlen(cfr->line) == 0)
		get_next_line(cfr);

	if (sscanf(cfr->line, "%s", var_name) != 1)
		goto EXIT;

	forward_pointer_in_the_line(cfr);
	status = 0;
EXIT:
	return status;	
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
