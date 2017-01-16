#include "nb/cfreader_bot.h"
#include "nb/geometric_bot.h"

typedef struct {
	void (*write_data)(FILE *fp, void *msh);
	int (*read_data)(nb_cfreader_t *cfr, void *msh);
} nbt_io;

static void nbt_io_init(nbt_io *io, const nb_mesh2D_t *mesh);
static void set_msh3trg_interface(nbt_io *io);
static void set_mshquad_interface(nbt_io *io);
static void set_mshpoly_interface(nbt_io *io);
static void set_mshpack_interface(nbt_io *io);
static void write_header(FILE *fp, const nb_mesh_t *mesh);
static int check_header(nb_cfreader_t *cfr);
static int read_type(nb_cfreader_t *cfr);
static void write_data(FILE *fp, const nb_mesh_t *mesh);

int nb_mesh2D_save_nbt(const nb_mesh2D_t *mesh,  const char *name)
{
	int status = 1;
	FILE *fp = fopen(name, "w");
	if (NULL == fp)
		goto EXIT;

	nbt_io io;
	nbt_io_init(&io, mesh);

	write_header(fp, mesh);
	io.write_data(fp, mesh->msh);

	status = 0;
CLOSE_FILE:
	fclose(fp);
EXIT:
	return status;
}

static void nbt_io_init(nbt_io *io, const nb_mesh2D_t *mesh)
{
	nb_mesh2D_type t = nb_mesh2D_get_type(mesh);
	const char *out;
	switch (t) {
	case NB_TRIAN:
		set_msh3trg_interface(io);
		break;
	case NB_QUAD:
		set_mshquad_interface(io);
		break;
	case NB_POLY:
		set_mshpoly_interface(io);
		break;
	case NB_DISK:
		set_mshpack_interface(io);
		break;
	}
}

static void set_msh3trg_interface(nbt_io *io)
{
	io->write_data = nb_msh3trg_write_data_nbt;
	io->read_data = nb_msh3trg_read_data_nbt;
}

static void set_mshquad_interface(nbt_io *io)
{
	io->write_data = nb_mshquad_write_data_nbt;
	io->read_data = nb_mshquad_read_data_nbt;
}

static void set_mshpoly_interface(nbt_io *io)
{
	io->write_data = nb_mshpoly_write_data_nbt;
	io->read_data = nb_mshpoly_read_data_nbt;
}

static void set_mshpack_interface(nbt_io *io)
{
	io->write_data = nb_mshpack_write_data_nbt;
	io->read_data = nb_mshpack_read_data_nbt;
}

static void write_header(FILE *fp, const nb_mesh_t *mesh)
{
	fprintf(fp, "[Numerical Bots File Format v1.0]\n");
	fprintf(fp, "Class = nb_mesh2D_t\n");
	const char *type = nb_mesh2D_get_type_string(mesh, type);
	fprintf(fp, "Type = %s\n", type);	
}

int nb_mesh2D_read_type_nbt(const char *name, nb_mesh2D_type *type)
{
	nb_cfreader_t *cfr = nb_cfreader_create();
	nb_cfreader_add_line_comment_token(cfr, "#");
	nb_cfreader_add_assignment_token(cfr, "=");
	int status = nb_cfreader_open_file(cfr, name);
	if (status != 0)
		goto EXIT;

	status = check_header(cfr);
	if (0 != status)
		goto CLOSE_FILE;

	status = read_type(cfr, type);
	if (0 != status)
		goto CLOSE_FILE;

	status = 0;
CLOSE_FILE:
	nb_cfreader_close_file(cfr);
EXIT:
	nb_cfreader_destroy(cfr);
	return status;
}

static int check_header(nb_cfreader_t *cfr)
{
	int status = nb_cfreader_check_line(cfr, "[Numerical Bots File Format v1.0]");
	if (0 != status)
		goto EXIT;

	status = nb_cfreader_check_line(cfr, "Class = nb_mesh2D_t");
EXIT:
	return status;
}

static int read_type(nb_cfreader_t *cfr)
{
	
EXIT:
	return status;
}

int nb_mesh2D_read_nbt(nb_mesh2D_t *mesh, const char *name)
{
	return 1;
}
