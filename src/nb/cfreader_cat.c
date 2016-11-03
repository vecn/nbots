/******************************************************************************
 *   CFreader CAT: Util to read customized plain text files                   *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "nb/memory_bot.h"

#include "nb/cfreader_cat.h"

#define LINE_LENGTH 256
#define TEMPORAL_TOKEN_LENGTH 10

struct nb_cfreader_s {
	/* OPPORTUNITY: Support several tokens for comments, multiline comments,
	 * open_close_string (of several chars, e.g. open/close string with << >>),
	 * etc
	 */
	char line_comment_token[TEMPORAL_TOKEN_LENGTH];
	char open_string_token;
	char close_string_token;
	FILE *fp;
	char buffer[LINE_LENGTH];
	char* line;
	uint32_t line_counter;
};

nb_cfreader_t* nb_cfreader_create(const char* filename,
				    const char* line_comment_token)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
		return NULL;
	nb_cfreader_t *cfreader = nb_allocate_mem(sizeof(*cfreader));
	sprintf(cfreader->line_comment_token, "%s", line_comment_token);
	cfreader->open_string_token  = '\"';
	cfreader->close_string_token = '\"';
	cfreader->fp = fp;
	cfreader->buffer[0] = '\0';
	cfreader->line = cfreader->buffer;
	cfreader->line_counter = 0;
	return cfreader;
}

char nb_cfreader_next_line(nb_cfreader_t* cfreader)
{
	while (fgets(cfreader->buffer, LINE_LENGTH, cfreader->fp) != NULL) {
		cfreader->line_counter += 1;
		/* Set end of line */
		char* pch = strchr(cfreader->buffer, '\n');
		if (pch != NULL)
			pch[0] = '\0';
		/* Handle line comments */
		pch = strstr(cfreader->buffer, cfreader->line_comment_token);
		if (pch != NULL)
			pch[0] = '\0';  /* Set end of line */
		/* Check if the line is empty */
		cfreader->line = cfreader->buffer;
		while (cfreader->line[0] == ' ' || cfreader->line[0] == '\t')
			cfreader->line = &(cfreader->line[1]);

		if (strlen(cfreader->line) == 0)
			continue;

		/* The line contains information */
		return 0;
	}
	return 1;
}

char nb_cfreader_read_int(nb_cfreader_t* cfreader, int* val)
{
	if (strlen(cfreader->line) == 0)
		nb_cfreader_next_line(cfreader);
	/* Read value */
	if (sscanf(cfreader->line, "%i", val) != 1)
		return 1;
	/* Forward char pointer in the line */
	while (cfreader->line[0] != ' ' &&
	       cfreader->line[0] != '\t' &&
	       cfreader->line[0] != '\0')
		cfreader->line = &(cfreader->line[1]);
	while (cfreader->line[0] == ' ' ||
	       cfreader->line[0] == '\t')
		cfreader->line = &(cfreader->line[1]);
	return 0;
}

char nb_cfreader_read_uint(nb_cfreader_t *cfreader, unsigned int *val)
{
	if (strlen(cfreader->line) == 0)
		nb_cfreader_next_line(cfreader);
	/* Read value */
	if (sscanf(cfreader->line, "%u", val) != 1)
		return 1;
	/* Forward char pointer in the line */
	while (cfreader->line[0] != ' ' &&
	       cfreader->line[0] != '\t' &&
	       cfreader->line[0] != '\0')
		cfreader->line = &(cfreader->line[1]);
	while (cfreader->line[0] == ' ' ||
	       cfreader->line[0] == '\t')
		cfreader->line = &(cfreader->line[1]);
	return 0;
}

char nb_cfreader_read_float(nb_cfreader_t* cfreader, float* val)
{
	if (strlen(cfreader->line) == 0)
		nb_cfreader_next_line(cfreader);
	/* Read value */
	if (sscanf(cfreader->line, "%f", val) != 1)
		return 1;
	/* Forward char pointer in the line */
	while (cfreader->line[0] != ' ' &&
	       cfreader->line[0] != '\t' &&
	       cfreader->line[0] != '\0')
		cfreader->line = &(cfreader->line[1]);
	while (cfreader->line[0] == ' ' ||
	       cfreader->line[0] == '\t')
		cfreader->line = &(cfreader->line[1]);
	return 0;
}
char nb_cfreader_read_double(nb_cfreader_t* cfreader, double* val)
{
	if (strlen(cfreader->line) == 0)
		nb_cfreader_next_line(cfreader);
	/* Read value */
	if (sscanf(cfreader->line, "%lf", val) != 1)
		return 1;
	/* Forward char pointer in the line */
	while (cfreader->line[0] != ' ' &&
	       cfreader->line[0] != '\t' &&
	       cfreader->line[0] != '\0')
		cfreader->line = &(cfreader->line[1]);
	while (cfreader->line[0] == ' ' ||
	       cfreader->line[0] == '\t')
		cfreader->line = &(cfreader->line[1]);
	return 0;
}

char nb_cfreader_read_bool(nb_cfreader_t* cfreader, bool* val)
{
	if(strlen(cfreader->line) == 0)
		nb_cfreader_next_line(cfreader);

	char* ptr = cfreader->line;
	uint32_t char_counter = 0;
	while (ptr[0] != ' ' && ptr[0] != '\t' && ptr[0] != '\0') {
		char_counter ++;
		ptr = &(ptr[1]);
	}
	char* value = nb_allocate_mem(char_counter+1);
	memcpy(value, cfreader->line, char_counter);
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

	/* Forward char pointer in the line */
	while (cfreader->line[0] != ' ' &&
	       cfreader->line[0] != '\t' &&
	       cfreader->line[0] != '\0')
		cfreader->line = &(cfreader->line[1]);
	while (cfreader->line[0] == ' ' || cfreader->line[0] == '\t')
		cfreader->line = &(cfreader->line[1]);

	/* Show warning of wrong formating */
	if (!is_a_correct_value) {
		printf("WARNING CFReader: ");
		printf("The boolean at line %i has an invalid value: \"%s\").\n",
		       cfreader->line_counter, value);
	}
	/* Free memory */
	nb_free_mem(value);

	/* Successful exit */
	return 0;
}

char* nb_cfreader_read_and_allocate_string(nb_cfreader_t* cfreader)
{
	if (strlen(cfreader->line) == 0)
		nb_cfreader_next_line(cfreader);
	/* Allocate and read string */
	bool throw_open_warning = false;
	if (cfreader->line[0] != cfreader->open_string_token)
		throw_open_warning = true;
	else
		cfreader->line = &(cfreader->line[1]);

	char* ptr = cfreader->line;
	uint32_t char_counter = 0;
	bool throw_close_warning = false;
	while (ptr[0] != cfreader->close_string_token) {
		if (ptr[0] == '\0') {
			throw_close_warning = true;
			break;
		}
		char_counter ++;
		ptr = &(ptr[1]);
	}

	char* string = NULL;
	if (char_counter > 0) {
		string = nb_allocate_mem(char_counter+1);
		memcpy(string, cfreader->line, char_counter);
		string[char_counter] = '\0';
		cfreader->line = &(cfreader->line[char_counter]);
		/* Throw warnings */
		if (throw_open_warning && throw_close_warning) {
			printf("WARNING CFReader: ");
			printf("Missing open/close tokens '%c/%c' of string '%s' at line %i.\n",
			       cfreader->open_string_token,
			       cfreader->close_string_token,
			       string,
			       cfreader->line_counter);
		} else if (throw_open_warning) {
			printf("WARNING CFReader: ");
			printf("Missing open token '%c' of string '%s' at line %i.\n",
			       cfreader->open_string_token,
			       string,
			       cfreader->line_counter);
		} else if (throw_close_warning) {
			printf("WARNING CFReader: ");
			printf("Missing close token '%c' of string '%s' at line %i.\n",
			       cfreader->close_string_token,
			       string,
			       cfreader->line_counter);
		}
	}
	if (ptr[0] == cfreader->close_string_token)
		cfreader->line = &(cfreader->line[1]);

	/* Forward char pointer in the line */
	while (cfreader->line[0] == ' ' || cfreader->line[0] == '\t')
		cfreader->line = &(cfreader->line[1]);

	/*  Return string */
	return string;
}

void nb_cfreader_destroy(nb_cfreader_t* cfreader)
{
	fclose(cfreader->fp);
	nb_free_mem(cfreader);
}
