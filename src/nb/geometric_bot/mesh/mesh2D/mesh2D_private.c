#include <stdlib.h>
#include <stdint.h>

#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/mesh2D/info.h"

#include "mesh2D_struct.h"
#include "mesh2D_private.h"

#include "elements2D/msh3trg_private.h"
#include "elements2D/mshquad_private.h"
#include "elements2D/mshpoly_private.h"
#include "elements2D/mshpack_private.h"

static void set_msh3trg_interface(nb_mesh2D_private_i *priv);
static void set_mshquad_interface(nb_mesh2D_private_i *priv);
static void set_mshpoly_interface(nb_mesh2D_private_i *priv);
static void set_mshpack_interface(nb_mesh2D_private_i *priv);

void nb_mesh2D_init_private_interface(nb_mesh2D_private_i *priv,
					 const nb_mesh2D_t *part)
{
	switch (part->type) {
	case NB_TRIAN:
		set_msh3trg_interface(priv);
		break;
	case NB_QUAD:
		set_mshquad_interface(priv);
		break;
	case NB_POLY:
		set_mshpoly_interface(priv);
		break;
	case NB_DISK:
		set_mshpack_interface(priv);
		break;
	default:
		set_msh3trg_interface(priv);
		break;
	}	
}

static void set_msh3trg_interface(nb_mesh2D_private_i *priv)
{
	priv->node_move_x = nb_msh3trg_node_move_x;
	priv->node_move_y = nb_msh3trg_node_move_y;
	priv->elem_move_x = nb_msh3trg_elem_move_x;
	priv->elem_move_y = nb_msh3trg_elem_move_y;
}

static void set_mshquad_interface(nb_mesh2D_private_i *priv)
{
	priv->node_move_x = nb_mshquad_node_move_x;
	priv->node_move_y = nb_mshquad_node_move_y;
	priv->elem_move_x = nb_mshquad_elem_move_x;
	priv->elem_move_y = nb_mshquad_elem_move_y;
}

static void set_mshpoly_interface(nb_mesh2D_private_i *priv)
{
	priv->node_move_x = nb_mshpoly_node_move_x;
	priv->node_move_y = nb_mshpoly_node_move_y;
	priv->elem_move_x = nb_mshpoly_elem_move_x;
	priv->elem_move_y = nb_mshpoly_elem_move_y;
}

static void set_mshpack_interface(nb_mesh2D_private_i *priv)
{
	priv->node_move_x = nb_mshpack_node_move_x;
	priv->node_move_y = nb_mshpack_node_move_y;
	priv->elem_move_x = nb_mshpack_elem_move_x;
	priv->elem_move_y = nb_mshpack_elem_move_y;
}
