#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "nb/memory_bot.h"
#include "nb/io_bot.h"
#include "nb/geometric_bot.h"

#include "mshquad_struct.h"
#include "mshquad_file_format_nbt.h"

void nb_mshquad_write_data_nbt(FILE *fp, const void *msh)
{
  ;
}


int nb_mshquad_read_data_nbt(nb_cfreader_t *cfr, void *msh_ptr)
{
  return 1;
}
