#ifndef __NB_IO_BOT_CFREADER_NBT_FILE_FORMAT_H__
#define __NB_IO_BOT_CFREADER_NBT_FILE_FORMAT_H__

#define NB_NBT_FILE_FORMAT_HEADER "[Numerical Bots File Format v1.0]"

void nb_cfreader_load_nbt_format(nb_cfreader_t *cfr);
int nb_cfreader_nbt_check_header(nb_cfreader_t *cfr, char *class);

#endif
