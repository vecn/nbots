#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "vcn/generic_dst.h"
#include "vcn/ugrid_dst.h"
#include "vcn/sparse_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/geometric_bot-cairo.h"
#include "geometric_bot-tests.c"

static double get_diff_seconds(struct timeval start, struct timeval end);

int main(int argc, char* argv[]){
  if(argc < 3){
    printf("Benchmark.Libre_Mesh must receive:\n");
    printf("   > Directory with the inputs\n");
    printf("   > Directory to save the outputs\n");
    printf("   > [OPTIONAL: Id of the test to perform]\n\n");
    return 0;
  }

  uint test_id = 0;
  if(argc > 3)
    test_id = atoi(argv[3]);
  
  /* Run benchmark */
  char name[100];
  struct timeval start, end;
  double time_cum = 0; /* Default initialization (Must be init) */

  for (register uint i = 0; i < N_tests; i++) {
    if (i != test_id - 1 && test_id != 0)
      continue;

    printf("Computing next test...");
    fflush(stdout);

    /* Run test */
    gettimeofday(&start, NULL);
    vcn_mesh_t*  mesh = test[i](name, argv[1]);
    gettimeofday(&end, NULL);

    if (mesh == NULL) {
      printf("\r%i: %s fails           \n", i+1, name);
      continue;
    }

    /* Show details */
    printf("\r%i: %s              \n", i+1, name);
    printf("\t Computing time: %lf\n\
\t Vertices: %i \n\
\t Edges: %i \n\
\t Triangles: %i \n", 
	   get_diff_seconds(start, end), 
	   vcn_mesh_get_number_of_vertices(mesh),
	   vcn_mesh_get_number_of_segments(mesh),
	   vcn_mesh_get_number_of_triangles(mesh));

    time_cum += get_diff_seconds(start, end);

    /* Define colors */
      double color_black[3] = {0.0, 0.0, 0.0};
      double color_white[4] = {1.0, 1.0, 1.0, 1.0};

    /* Visualize Delaunay Mesh */
    vcn_msh3trg_t* msh3trg =
      vcn_mesh_get_msh3trg(mesh, false, true, true, true, true, NULL);
    sprintf(name, "%s/test%i_msh3trg.png", argv[2], i+1);
    vcn_msh3trg_save_png(msh3trg, name, 1000, 800, 
			 color_white, color_white, 
			 color_black, color_black, 0.5, 1.0);

    /* Visualize partitions */
    if (true) {
      uint N_part = 16;
      vcn_graph_t* graph = vcn_msh3trg_create_elem_graph(msh3trg);
      uint* part = vcn_graph_partition_sb(graph, N_part);
      sprintf(name, "%s/test%i_msh3trg_partition.png", argv[2], i+1);
      double color_line[3] = {0.0, 0.0, 0.0};
      vcn_msh3trg_partition_save_png(msh3trg, name, 1000, 800, 
				     N_part, part, N_part, 1.0,
				     NULL, 0.7, color_line, color_line,
				     0.5, 5.0, false);
      free(part);
    }
    vcn_msh3trg_destroy(msh3trg);

    /* Visualize Voronoi Diagram */
    if (false) {
      vcn_mshpoly_t* mshpoly =
	vcn_mesh_get_mshpoly(mesh, false, false, 0, NULL, NULL);
      sprintf(name, "%s/test%i_mshpoly.png", argv[2], i+1);
      vcn_mshpoly_save_png(mshpoly, name, false, 1000, 800);
      vcn_mshpoly_destroy(mshpoly);
      
      /* Visualize Voronoi Diagram after Lloyd */
      if (true) {
	mshpoly = vcn_mesh_get_mshpoly(mesh, false, true, 100, NULL, NULL);
	sprintf(name, "%s/test%i_mshpoly_central.png", argv[2], i+1);
	vcn_mshpoly_save_png(mshpoly, name, false, 1000, 800);
	vcn_mshpoly_destroy(mshpoly);
      }
    }
    
    /* Visualize Sphere Packing */
    if (false) {
      vcn_mshpack_t* mshpack = vcn_mesh_get_mshpack(mesh, false, 100,
						    0.1, 0.1, NULL);
      sprintf(name, "%s/test%i_mshpack.png", argv[2], i+1);
      vcn_mshpack_save_png(mshpack, name, 1000, 800);
      vcn_mshpack_destroy(mshpack);
    }

    /* Destroy mesh */
    vcn_mesh_destroy(mesh);
  }
  /* Finish */
  printf("\nTotal time: %lf\n", time_cum);

  /* Successful finish */
  return 0;
}

static double get_diff_seconds(struct timeval start, struct timeval end){
  return (double)((end.tv_sec*1000000 + end.tv_usec) -
		  (start.tv_sec*1000000 + start.tv_usec)) / 1e6;
}
