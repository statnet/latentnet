#ifndef GEODIST_H
#define GEODIST_H

#include "edgeTree.h"

void node_geodesics (int *edgelist, int *nnodes, int *nodelist,
                     int *nedges, int *nodecolor, int *dist, 
                     int *Q, int *source);

void geodesic_matrix (int *edgelist, int *nnodes,
		      int *nodelist, int *nedges, 
		      int *nodecolor, int *distmat, int *Q);

#endif
