#include "nb/io_bot.h"

void nb_cfreader_load_nbt_format(nb_cfreader_t *cfr)
{
	nb_cfreader_add_line_comment_token(cfr, "#");
	nb_cfreader_add_line_comment_token(cfr, "//");
	nb_cfreader_add_assignment_token(cfr, "=");
	nb_cfreader_add_assignment_token(cfr, "<-");
}

int nb_cfreader_nbt_check_header(nb_cfreader_t *cfr, char *class)
{
	int status = nb_cfreader_check_line(cfr, NB_NBT_FILE_FORMAT_HEADER);
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_read_var_token(cfr, "Class", class);
EXIT:
	return status;

}
