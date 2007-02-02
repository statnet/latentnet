#ifndef EDGETREE_H
#define EDGETREE_H

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MAXEDGES 100000
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))

typedef unsigned long Vertex;
typedef unsigned long Edge;

/*  next_inedge and next_outedge are continually updated to give
    the smallest index of an edge object not being used.
*/
Edge next_inedge;
Edge next_outedge;

Edge next_minedge;
Edge next_moutedge;

Edge next_inedge0;
Edge next_outedge0;

/*  outdegree[] and indegree[] are continually updated to give
    the appropriate degree values for each vertex.  These should
    point to vectors of Vertex of length nvertices.  n_nodes is
    simply the number of vertices in the graph.  directed_flag is
    1 or 0, depending on whether or not the graph is directed.
*/
Vertex *outdegree;
Vertex *indegree;
Vertex n_nodes;
Vertex n_edges;
int directed_flag;

Vertex designempty;
Vertex *moutdegree;
Vertex *mindegree;
Vertex n_medges;

Vertex *outdegree0;
Vertex *indegree0;
Vertex n_edges0;

/*  Edgestruct is a binary tree structure, which is how the edgelists 
    are stored.  The root of the tree for vertex i will be inedges[i]
    or outedges[i].  inedges[0] and outedges[0] are unused, since the
    index 0 will indicate no link.  Indices are long unsigned integers,
    which means graphs can contain 2^32-1= 4294967295 edges (enough to
    hold all edges in a 92682-vertex undirected graph or a 65536-vertex
    directed graph, assuming no multiple edges or loops), though for this
    MAXEDGES must be adjusted accordingly.
*/
typedef struct estruct {
  Vertex value;      /*  the other end of the edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
} Edgestruct;

/*  inedges and outedges are arrays of Edgestruct that are used to 
    store all of the incoming and outgoing edges, respectively.
*/
Edgestruct inedges[MAXEDGES], outedges[MAXEDGES];
Edgestruct minedges[MAXEDGES], moutedges[MAXEDGES];
Edgestruct inedges0[MAXEDGES], outedges0[MAXEDGES];

/* The Gptr is a structure containing pointers to all essential elements
   of a given graph; it is intended to be used to allow the various
   edgeTree routines to alter the same graph, even when being called
   in a foreign namespace (in which global variables would be out of
   scope.
*/
typedef struct Gptrtype {
  Edgestruct *inedges;
  Edgestruct *outedges;
  int *directed_flag;
  Vertex *n_nodes;
  Vertex *n_edges;
  Edge *next_inedge;
  Edge *next_outedge;
  Vertex *indegree;
  Vertex *outdegree;

/*  The Design Variables */

  Vertex *designempty;
  Edgestruct *minedges;
  Edgestruct *moutedges;
  Vertex *n_medges;
  Edge *next_minedge;
  Edge *next_moutedge;
  Vertex *mindegree;
  Vertex *moutdegree;
} Gptr;

Gptr GraphInitialize(double *heads, double *tails, Edge nedges);
Gptr DesignInitialize(double *heads, double *tails, Edge nedges);
Gptr GraphInitialize0(double *heads, double *tails, Edge nedges);
int DesignMissing (Vertex a, Vertex b, Gptr mg);
void GraphDestroy();
void GraphDestroy0();
void DesignDestroy();

Edge EdgetreeSearch (Vertex a, Vertex b, Edgestruct *edges);
Edge EdgetreeSuccessor (Edgestruct *edges, Edge x);
Edge EdgetreeMinimum (Edgestruct *edges, Edge x);
int EdgeToValue(Vertex head, Vertex tail);
int ToggleEdge (Vertex head, Vertex tail, Gptr g);
int AddEdgeToTrees(Vertex head, Vertex tail, Gptr g);
void AddHalfedgeToTree (Vertex a, Vertex b, Edgestruct *edges, 
		   Edge *next_edge);
int DeleteEdgeFromTrees(Vertex head, Vertex tail, Gptr g);
int DeleteHalfedgeFromTree(Vertex a, Vertex b, Edgestruct *edges,
		     Edge *next_edge);
void printedge(Edge e, Edgestruct *edges);
void InOrderTreeWalk(Edgestruct *edges, Edge x);

#endif
