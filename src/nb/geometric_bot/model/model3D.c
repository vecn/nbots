/* Por implementar, responsable: Jorge Lopez */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#include "nb/memory_bot/allocate_mem.h"
#include "nb/geometric_bot/model/model3D_struct.h"
#include "nb/geometric_bot/model/model3D.h"

uint64_t nb_model3D_get_memsize( uint32_t n_triangles )
{
        uint64_t memsize = 0;
        memsize += sizeof( nb_model3D_t );
        memsize += ( 9 * n_triangles * sizeof( double ) );
        memsize += ( 3 * n_triangles * sizeof( double ) );
        memsize += ( 3 * n_triangles * sizeof( uint32_t ) );
        if( memsize == 0 ){
                printf(" \n");
                printf("   /************************************************************/\n" );
                printf("   /*No memory calculated for boundary mesh,finishing algorithm*/\n" );
                printf("   /************************************************************/\n\n" );
                assert(memsize);
        }
        return memsize;
}

void nb_model3D_init(void *model_ptr)
{
  ;
}

void nb_model3D_copy(void *model_ptr, const void *src_model_ptr)
{
        nb_model3D_t *src_model,*dst_model;
        src_model = src_model_ptr;
        dst_model = model_ptr;
        uint32_t n_triangles = nb_model3D_get_N_face( src_model );
        uint32_t n_vtx = nb_model3D_get_N_vtx( src_model );
        nb_model3D_set_N_vtx ( n_vtx , dst_model );
        nb_model3D_set_N_face( n_triangles , dst_model );
        for( int i_vertex = 0 ; i_vertex < ( n_triangles * 9 ) ; i_vertex++ ){
                double vtx = nb_model3D_get_vtx( i_vertex , src_model );
                nb_model3D_set_vtx( i_vertex , vtx , dst_model );
        }
        for( int i_normal = 0 ; i_normal < ( n_triangles * 3 ) ; i_normal++ ){
                double nf = nb_model3D_get_nf( i_normal , src_model );
                nb_model3D_set_nf( i_normal , nf , dst_model );
        }
        for( int i_adj = 0 ; i_adj < ( n_triangles * 3 ) ; i_adj++ ){
                uint32_t adj = nb_model3D_get_adj( i_adj , src_model );
                nb_model3D_set_adj( i_adj , adj , dst_model );
        }

}

void nb_model3D_finish(void *model_ptr)
{
	;
}

void nb_model3D_clear(void *model_ptr)
{
        free(model_ptr);
        model_ptr = NULL;
}

void nb_model3D_load(void *model_ptr, const char* filename, uint32_t n_triangles)
{
        nb_model3D_t *model = model_ptr;
        nb_model3D_load_model_from_ASCII_STL( filename , model , n_triangles);
}

void nb_model3D_load_box(void *model_ptr,
			 double x_min, double y_min, double z_min,
			 double x_max, double y_max, double z_max)
{
	;
}

int nb_model3D_save(const void *model, const char* filename)
{
  FILE *fp;
  fp=fopen(filename,"w");
  nb_model3D_t *src_model;
  src_model = model;
  fprintf(fp,"solid name\n");
  uint32_t n_faces = nb_model3D_get_N_face( src_model );
  for(int i_face = 0 ; i_face < n_faces ; i_face++ ){
    nb_model3D_print_face_in_file( fp , i_face , src_model );
  }
  fprintf(fp,"endsolid name\n");
  fclose(fp);
  return 1;
}

void nb_model3D_get_enveloping_box(const void *const model, double box[6])
{

  nb_model3D_t *src_model;
  src_model = model;
  box[ 0 ] = 1e20;    box[ 1 ] = -1e20;
  box[ 2 ] = 1e20;    box[ 3 ] = -1e20;
  box[ 4 ] = 1e20;    box[ 5 ] = -1e20;
  uint32_t n_faces = nb_model3D_get_N_face( src_model );
  uint32_t n_coords = n_faces * 3;
  for(  int i_coord = 0  ;  i_coord < n_coords  ;  i_coord++  ){
        double x_coord = nb_model3D_get_vtx( ( i_coord * 3 ), src_model );
        double y_coord = nb_model3D_get_vtx( ( i_coord * 3 ) + 1, src_model );
        double z_coord = nb_model3D_get_vtx( ( i_coord * 3 ) + 2, src_model );
        if( x_coord < box[ 0 ] ){
          box[ 0 ] = x_coord;
        }
        if( x_coord > box[ 1 ] ){
          box[ 1 ] = x_coord;
        }
        if( y_coord < box[ 2 ] ){
          box[ 2 ] = y_coord;
        }
        if( y_coord > box[ 3 ] ){
          box[ 3 ] = y_coord;
        }
        if( z_coord < box[ 4 ] ){
          box[ 4 ] = z_coord;
        }
        if( z_coord > box[ 5 ] ){
          box[ 5 ] = z_coord;
        }
  }
}

uint32_t nb_model3D_get_N_vtx(const void *const model)
{
        nb_model3D_t *model_t = model;
        return model_t->N_vtx;
}

uint32_t nb_model3D_get_N_edges(const void *const model)
{
	;
}

uint32_t nb_model3D_get_N_face(const void *const model)
{
        nb_model3D_t *model_t = model;
        return model_t->N_face;
}


uint32_t nb_model3D_get_number_of_triangles_in_ASCII_STL( const char *name ){
        FILE *fp = NULL;
        char header_name[ 100 ];
        int cont = 0;
        fp = fopen( name , "r" );
        assert( fp );
        nb_model3D_read_header_ASCII_STL( fp , header_name );
        while( 1 ){
                char aux_read[ 100 ];
                fscanf( fp , "%s" , aux_read );
                if( !strcmp( aux_read , "endsolid\0" ) ){
                        break;
                }
                while( 1 ){
                        fscanf( fp , "%s" , aux_read );
                        if( !strcmp( aux_read , "endfacet\0" ) ){
                                break;
                        }
                }
                cont++;
        }
        fclose(fp);
        if( cont == 0 ){
                printf(" \n");
                printf("   /****************************************************/\n" );
                printf("   /*No triangles on boundary mesh, finishing algorithm*/\n" );
                printf("   /****************************************************/\n\n" );
                assert(cont);
        }
        return cont;
}

void nb_model3D_read_header_ASCII_STL( FILE *fp , char *header ){
        char read[ 100 ];
        fscanf( fp , "%s" , read );
        fscanf( fp , "%s" , read );
        strcpy( header , read );
}


nb_model3D_t* nb_model3D_assign_memory_on_struct_for_ASCII_STL(uint64_t memsize , 
                                                               uint32_t n_triangles){
        nb_model3D_t *model = NULL;
        char *ptr = NULL;
        ptr = nb_allocate_mem( memsize );
        if( ptr == NULL){
                printf(" \n");
                printf("   /***********************************************************/\n" );
                printf("   /*No memory allocated for boundary mesh,finishing algorithm*/\n" );
                printf("   /***********************************************************/\n\n" );
                assert( ptr );
        }
        model = ptr;
        model->vtx = &ptr[ 40 ];
        model->nf  = &ptr[ 40 + ( n_triangles * 9 * 8 ) ];
        model->adj = &ptr[ 40 + ( n_triangles * 9 * 8 ) + ( n_triangles * 3 * 8 ) ];
        return model;
}

void nb_model3D_load_model_from_ASCII_STL( const char *name , nb_model3D_t *boundary_t ,
                                           uint32_t n_triangles ){
        int i_vertex = 0;
        int i_normal = 0;
        int i_adj    = 0;
        double aux_double;
        FILE *fp = NULL;
        char header_name[ 100 ];
        fp = fopen( name , "r" );
        assert( fp );
        nb_model3D_read_header_ASCII_STL( fp , header_name );
        while( 1 ){
                char aux_read[ 100 ];
                fscanf( fp , "%s" , aux_read );
                if( !strcmp( aux_read , "endsolid\0" ) ){
                        break;
                }
                fscanf( fp , "%s" , aux_read );
                fscanf( fp , "%lf" , &aux_double );
                nb_model3D_set_nf( i_normal , aux_double , boundary_t ); i_normal++;
                fscanf( fp , "%lf" , &aux_double );
                nb_model3D_set_nf( i_normal , aux_double , boundary_t ); i_normal++;
                fscanf( fp , "%lf" , &aux_double );
                nb_model3D_set_nf( i_normal , aux_double , boundary_t ); i_normal++;
                fscanf( fp , "%s" , aux_read );    fscanf( fp , "%s" , aux_read );

                for(int j_vertex = 0 ; j_vertex<3;j_vertex++){
                        fscanf( fp , "%s" , aux_read );
                        fscanf( fp , "%lf" , &aux_double );
                        nb_model3D_set_vtx( i_vertex , aux_double , boundary_t ); i_vertex++;
                        fscanf( fp , "%lf" , &aux_double );
                        nb_model3D_set_vtx( i_vertex , aux_double , boundary_t ); i_vertex++;
                        fscanf( fp , "%lf" , &aux_double );
                        nb_model3D_set_vtx( i_vertex , aux_double , boundary_t ); i_vertex++;
                }
                fscanf( fp , "%s" , aux_read );
                fscanf( fp , "%s" , aux_read );

                nb_model3D_set_adj( i_adj , (uint32_t)i_adj , boundary_t ); i_adj++;
                nb_model3D_set_adj( i_adj , (uint32_t)i_adj , boundary_t ); i_adj++;
                nb_model3D_set_adj( i_adj , (uint32_t)i_adj , boundary_t ); i_adj++;
        }
        fclose(fp);
        if( i_vertex == 0  ||  i_normal == 0  ||  i_adj == 0 ){
                printf(" \n");
                printf("   /****************************************************/\n" );
                printf("   /*No triangles on boundary mesh, finishing algorithm*/\n" );
                printf("   /****************************************************/\n\n" );
                assert(i_vertex);
                assert(i_normal);
                assert(i_adj   );
        }
        nb_model3D_set_N_vtx((n_triangles * 3),boundary_t);
        nb_model3D_set_N_face(n_triangles,boundary_t);
}

void nb_model3D_set_nf(int position,double value, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        model->nf[ position ] = value;
}

void nb_model3D_set_vtx(int position,double value, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        model->vtx[ position ] = value;
}

void nb_model3D_set_adj(int position,uint32_t value, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        model->adj[ position ] = value;
}

void nb_model3D_set_N_vtx(uint32_t value, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        model->N_vtx = value;
}

void nb_model3D_set_N_face(uint32_t value, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        model->N_face = value;
}

double nb_model3D_get_nf(int position, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        return model->nf[ position ];
}

double nb_model3D_get_vtx(int position, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        return model->vtx[ position ];
}

uint32_t nb_model3D_get_adj(int position, void *boundary_t){
        nb_model3D_t *model = boundary_t;
        return model->adj[ position ];
}

void nb_model3D_print_face_in_file( FILE *fp , int i_face , const void *model ){
        nb_model3D_t *src_model;
        src_model = model;
        fprintf(fp,"\tfacet normal\t");
        for(int i_normal = 0  ;  i_normal < 3  ;  i_normal++  ){
          int position = ( ( i_face * 3 ) + i_normal );
          double normal = nb_model3D_get_nf( position , src_model );
          fprintf( fp , "%lf\t" , normal );
        }
        fprintf(fp,"\n\t\touter loop\n");
        for(  int i_vertex = 0  ;  i_vertex < 3  ;  i_vertex++  ){
          fprintf( fp , "\t\t\tvertex\t" );
          int position = ( i_face * 3) + i_vertex;
          uint32_t coord_index = nb_model3D_get_adj( position , src_model );
          for(  int i_coord=0  ;  i_coord < 3  ;  i_coord++  ){
            int coord_position = ( ( (int)coord_index ) * 3 ) + i_coord;
            double coord = nb_model3D_get_vtx( coord_position , src_model );
            fprintf( fp , "%lf\t" , coord );
          }
          fprintf(fp,"\n");
        }
        fprintf(fp,"\t\tendloop\n");
        fprintf(fp,"\tendfacet\n");
}










































