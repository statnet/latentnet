#include "edgeTree.h"


/*  ***************** */
/*  void GraphInitialize */
/*  */
/*  Initialize, construct binary tree version of graph.  Note */
/*  that zero Edgestruct should have all its values set to zero */
/*  ***************** */

Gptr GraphInitialize(double *heads, double *tails, Edge nedges) {
  Edge i, j;
  Vertex h, t;
  Gptr g;

  /* Initialize next_edge's, allocate space for degree vectors */
  next_inedge = next_outedge = (Edge)n_nodes+1;
  outdegree = (Vertex *) malloc(sizeof(Vertex) * (n_nodes+1));
  indegree  = (Vertex *) malloc(sizeof(Vertex) * (n_nodes+1));
  
  for (i=0; i<=n_nodes; i++) {
    inedges[i].value = outedges[i].value = 0;
    inedges[i].parent = outedges[i].parent = 0;
    inedges[i].left = outedges[i].left = 0;
    inedges[i].right = outedges[i].right = 0;
    outdegree[i] = indegree[i] = 0;
  }
  
  for (; i<MAXEDGES; i++)
    inedges[i].value = outedges[i].value = 0;

  /*Configure a Gptr*/
  g.outedges=outedges;
  g.inedges=inedges;
  g.indegree=indegree;
  g.outdegree=outdegree;
  g.n_nodes=&n_nodes;
  g.n_edges=&n_edges;
  g.next_inedge=&next_inedge;
  g.next_outedge=&next_outedge;
  g.directed_flag=&directed_flag;

  *g.n_edges=*g.n_edges - n_edges;
  for(i = nedges; i > 0; i--) {
    j = i * unif_rand();  /* shuffle edgelist to help bin. tree achieve best perf */
    h = (Vertex)heads[j];
    t = (Vertex)tails[j];
    heads[j] = heads[i-1];
    tails[j] = tails[i-1];
    heads[i-1] = h;
    tails[i-1] = t;
    if (!directed_flag && h > t) 
      AddEdgeToTrees(t,h,g); /* Undir edges always have head < tail */ 
    else 
      AddEdgeToTrees(h,t,g);
  }

#if 0
  Rprintf("Graph:\n");
  for (i = 0; i<50 ; i++)
    {
      Rprintf("%d  ",  inedges[i].value);
      if ((i+1)%5 == 0)
	Rprintf("\n");
    }
#endif

  return g;
}

/*  ***************** */
/*  void GraphDestroy */
/*  */
/*  ***************** */
void GraphDestroy() {
  int i;

  for (i = 0; i < MAXEDGES; i++) {
    inedges[i].value = outedges[i].value = 0;  
    inedges[i].parent = outedges[i].parent = 0;  
    inedges[i].left = outedges[i].left = 0;  
    inedges[i].right = outedges[i].right = 0;  
  }

  for (i=0; i<=n_nodes; i++)
    outdegree[i] = indegree[i] = 0;

  free(outdegree);
  free(indegree);  
}

/*  ***************** */
/*  Edge EdgetreeSearch */
/*  */
/*  Check to see if there's an Edgestruct with value b  */
/*  in the tree rooted at edges[a].  Return i such that  */
/*  edges[i] is that Edgestruct, or 0 if none. */
/*  ***************** */

Edge EdgetreeSearch (Vertex a, Vertex b, Edgestruct *edges) {
  Edgestruct *es;
  Edge e = a;
  Vertex v;
  /* this part is to try to sort out the below */
  es = edges + e;
  v = es->value;

  while(e != 0 && b != v)
    {
      e = (b<v)?  es->left : es->right;
      es = edges + e;
      v = es->value;
    }
  return e;

  /*  while (e!=0 && b!= (v= (es=edges+e)->value ) )
      e = (b<v)?  es->left : es->right; */
}

/*  ***************** */
/*  Edge EdgetreeSuccessor */
/*  */
/*  Return the index of the Edgestruct with the smallest value */
/*  greater than edges[z].value in the same edge tree, or 0 */
/*  if none.  This is not very useful directly, but it is  */
/*  called by the DeleteHalfedgeFromTree function. */
/*  ***************** */

Edge EdgetreeSuccessor (Edgestruct *edges, Edge x) {
  Edgestruct *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->right) != 0) 
    return EdgetreeMinimum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->right) 
    x=y;
  return y; 
}   

/*  ***************** */
/*  Edge EdgetreeMinimum */
/*  */
/*  Return the index of the Edgestruct with the */
/*  smallest value in the subtree rooted at x */
/*  ***************** */

Edge EdgetreeMinimum (Edgestruct *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}

/*  ***************** */
/*  Edge ToggleEdge */
/*  */
/*  Toggle an edge:  Set it to the opposite of its current */
/*  value.  Return 1 if edge added, 0 if deleted. */
/*  ***************** */

int ToggleEdge (Vertex head, Vertex tail, Gptr g) 
{
  if (!*(g.directed_flag) && head > tail) 
    {
      Vertex temp;
      temp = head; /*  Make sure head<tail always for undirected edges */
      head = tail;
      tail = temp;
    }

  if (AddEdgeToTrees(head,tail,g))
    return 1;
  else 
    return 1 - DeleteEdgeFromTrees(head,tail,g);
}

/*  ***************** */
/*  Edge AddEdgeToTrees */
/*  */
/*  Add an edge from head to tail after checking to see */
/*  if it's legal. Return 1 if edge added, 0 otherwise.   */
/*  ***************** */

int AddEdgeToTrees(Vertex head, Vertex tail, Gptr g){
  if (EdgetreeSearch(head, tail, g.outedges) == 0) {
    AddHalfedgeToTree(head, tail, g.outedges, g.next_outedge);
    AddHalfedgeToTree(tail, head, g.inedges, g.next_inedge);
    ++g.outdegree[head];
    ++g.indegree[tail];
    ++(*g.n_edges);
    return 1;
  }
  return 0;
}

/*  ***************** */
/*  Edge AddHalfedgeToTrees */
/*  */
/*  An "edge" should be added to both the list of outedges and the */
/*  list of inedges.  Thus, the AddEdgeToTrees routine calls AddHalfedgeToTree */
/*  twice, once in each direction. */
/*  ***************** */

void AddHalfedgeToTree (Vertex a, Vertex b, Edgestruct *edges, 
		   Edge *next_edge) 
{
  Edgestruct *eptr = edges+a, *new;
  Edge e;
  
  if (eptr->value==0) 
    { /* This is the first edge for this vertex. */
      eptr->value=b;
      return;
    }
  (new=edges+*next_edge)->value=b;  
  new->left=new->right=0;
  for (e=a; e!=0; e=(b < (eptr=edges+e)->value) ? eptr->left : eptr->right);
  new->parent=eptr-edges;
  if (b < eptr->value)
    eptr->left=*next_edge;
  else
    eptr->right=*next_edge;
  while (++*next_edge<MAXEDGES && edges[*next_edge].value!=0);
}

/*  ***************** */
/*  int DeleteEdgeFromTrees */
/*  */
/*  Find and delete the edge from head to tail.   */
/*  Return 1 if successful, 0 otherwise.   */
/*  ***************** */

int DeleteEdgeFromTrees(Vertex head, Vertex tail, Gptr g){
  if (DeleteHalfedgeFromTree(head, tail, g.outedges, g.next_outedge) &&
      DeleteHalfedgeFromTree(tail, head, g.inedges, g.next_inedge)) {
    --g.outdegree[head];
    --g.indegree[tail];
    --(*g.n_edges);
    return 1;
  }
  return 0;
}

/*  ***************** */
/*  int DeleteHalfedgeFromTree */
/*  */
/*  Delete the Edgestruct with value b from the tree rooted at edges[a]. */
/*  Return 0 if no such Edgestruct exists, 1 otherwise.  Also update the */
/*  value of *next_edge appropriately. */
/*  ***************** */

int DeleteHalfedgeFromTree(Vertex a, Vertex b, Edgestruct *edges,
		     Edge *next_edge){
  Edge x, z, root=(Edge)a;
  Edgestruct *xptr, *zptr, *ptr;

  if ((z=EdgetreeSearch(a, b, edges))==0)  /* z is the current Edgestruct. */
    return 0;
  /* First, determine which node to splice out; this is z */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    z=EdgetreeSuccessor(edges, z);
    zptr->value = (ptr=edges+z)->value;
    zptr=ptr;
  }
  /* Set x to the child of z (there is at most one). */
  if ((x=zptr->left) == 0)
    x = zptr->right;
  /* Splice out node z */
  if (z == root) {
    zptr->value = (xptr=edges+x)->value;
    if (x != 0) {
      if ((zptr->left=xptr->left) != 0)
	(edges+zptr->left)->parent = z;
      if ((zptr->right=xptr->right) != 0)
	(edges+zptr->right)->parent = z;
      zptr=edges+(z=x);
    }  else 
      return 1;
  } else {
    if (x != 0)
      (xptr=edges+x)->parent = zptr->parent;
    if (z==(ptr=(edges+zptr->parent))->left)
      ptr->left = x;
    else 
      ptr->right = x;
  }  
  /* Clear z node, update *next_edge if necessary. */
  zptr->value=0;
  if (z < *next_edge)
    *next_edge=z;
  return 1;
}

/*  ***************** */
/*  void printedge */
/*  */
/*  Diagnostic routine that prints out the contents */
/*  of the specified Edgestruct (used for debugging).   */
/*  ***************** */

void printedge(Edge e, Edgestruct *edges){
  Rprintf("Edge structure [%d]:\n",e);
  Rprintf("\t.value=%d\n",edges[e].value);
  Rprintf("\t.parent=%d\n",edges[e].parent);
  Rprintf("\t.left=%d\n",edges[e].left);
  Rprintf("\t.right=%d\n",edges[e].right);
}

/*  ***************** */
/*  void InOrderTreeWalk */
/*  */
/*  Diagnostic routine that prints the nodes in the tree rooted */
/*  at edges[x], in increasing order of their values. */
/*  ***************** */

void InOrderTreeWalk(Edgestruct *edges, Edge x) {
  if (x != 0) {
    InOrderTreeWalk(edges, (edges+x)->left);
    /*    printedge(x, edges); */
    Rprintf(" %d ",(edges+x)->value); 
    InOrderTreeWalk(edges, (edges+x)->right);
  }
}
Gptr DesignInitialize(double *heads, double *tails, Edge nedges) {
  Edge i, j;
  Vertex h, t;
  Gptr mg;

  /* Initialize next_edge's, allocate space for degree vectors */
  next_minedge = next_moutedge = (Edge)n_nodes+1;
  moutdegree = (Vertex *) malloc(sizeof(Vertex) * (n_nodes+1));
  mindegree  = (Vertex *) malloc(sizeof(Vertex) * (n_nodes+1));
  
  for (i=0; i<=n_nodes; i++) {
    minedges[i].value = moutedges[i].value = 0;
    minedges[i].parent = moutedges[i].parent = 0;
    minedges[i].left = moutedges[i].left = 0;
    minedges[i].right = moutedges[i].right = 0;
    moutdegree[i] = mindegree[i] = 0;
  }
  
  for (; i<MAXEDGES; i++)
    minedges[i].value = moutedges[i].value = 0;

  /*Configure a Gptr*/
  mg.designempty=&designempty;
  mg.outedges=moutedges;
  mg.inedges=minedges;
  mg.indegree=mindegree;
  mg.outdegree=moutdegree;
  mg.n_nodes=&n_nodes;
  mg.n_edges=&n_medges;
  mg.next_inedge=&next_minedge;
  mg.next_outedge=&next_moutedge;
  mg.directed_flag=&directed_flag;

  if(nedges > 0){
   mg.designempty = 0;
   for(i = nedges; i > 0; i--) {
    j = i * unif_rand();  /* shuffle edgelist to help bin. tree achieve best perf */
    h = (Vertex)heads[j];
    t = (Vertex)tails[j];
    heads[j] = heads[i-1];
    tails[j] = tails[i-1];
    heads[i-1] = h;
    tails[i-1] = t;
    if (!directed_flag && h > t) 
      AddEdgeToTrees(t,h,mg); /* Undir edges always have head < tail */ 
    else 
      AddEdgeToTrees(h,t,mg);
   }
  }else{
   mg.designempty = 1;
  }

#if 0
  Rprintf("Graph:\n");
  for (i = 0; i<50 ; i++)
    {
      Rprintf("%d  ",  minedges[i].value);
      if ((i+1)%5 == 0)
	Rprintf("\n");
    }
#endif

  return mg;
}

/*  ***************** */
/*  void DesignDestroy */
/*  */
/*  ***************** */
void DesignDestroy() {
  int i;

  for (i = 0; i < MAXEDGES; i++) {
    minedges[i].value = moutedges[i].value = 0;  
    minedges[i].parent = moutedges[i].parent = 0;  
    minedges[i].left = moutedges[i].left = 0;  
    minedges[i].right = moutedges[i].right = 0;  
  }

  for (i=0; i<=n_nodes; i++)
    moutdegree[i] = mindegree[i] = 0;

  free(moutdegree);
  free(mindegree);  
}

/*  ***************** */
/*  Edge EdgetreeSearch */
/*  */
/*  Check to see if there's an Edgestruct with value b  */
/*  in the tree rooted at edges[a].  Return i such that  */
/*  edges[i] is that Edgestruct, or 0 if none. */
/*  ***************** */

int DesignMissing (Vertex a, Vertex b, Gptr mg) {
/*   if(mg.designempty){return(1);} */
  int miss;
  miss = EdgetreeSearch(a,b,mg.outedges);
  if(*(mg.directed_flag)){
    miss += EdgetreeSearch(a,b,mg.inedges);
  }
  return(miss);
}
Gptr GraphInitialize0(double *heads, double *tails, Edge nedges) {
  Edge i, j;
  Vertex h, t;
  Gptr g0;

  /* Initialize next_edge's, allocate space for degree vectors */
  next_inedge0 = next_outedge0 = (Edge)n_nodes+1;
  outdegree0 = (Vertex *) malloc(sizeof(Vertex) * (n_nodes+1));
  indegree0  = (Vertex *) malloc(sizeof(Vertex) * (n_nodes+1));
  
  for (i=0; i<=n_nodes; i++) {
    inedges0[i].value = outedges0[i].value = 0;
    inedges0[i].parent = outedges0[i].parent = 0;
    inedges0[i].left = outedges0[i].left = 0;
    inedges0[i].right = outedges0[i].right = 0;
    outdegree0[i] = indegree0[i] = 0;
  }
  
  for (; i<MAXEDGES; i++)
    inedges0[i].value = outedges0[i].value = 0;

  /*Configure a Gptr*/
  g0.outedges=outedges0;
  g0.inedges=inedges0;
  g0.indegree=indegree0;
  g0.outdegree=outdegree0;
  g0.n_nodes=&n_nodes;
  g0.n_edges=&n_edges0;
  g0.next_inedge=&next_inedge0;
  g0.next_outedge=&next_outedge0;
  g0.directed_flag=&directed_flag;

  *g0.n_edges=*g0.n_edges - n_edges0;
  for(i = nedges; i > 0; i--) {
    j = i * unif_rand();  /* shuffle edgelist to help bin. tree achieve best perf */
    h = (Vertex)heads[j];
    t = (Vertex)tails[j];
    heads[j] = heads[i-1];
    tails[j] = tails[i-1];
    heads[i-1] = h;
    tails[i-1] = t;
    if (!directed_flag && h > t) 
      AddEdgeToTrees(t,h,g0); /* Undir edges always have head < tail */ 
    else 
      AddEdgeToTrees(h,t,g0);
  }

#if 0
  Rprintf("Graph:\n");
  for (i = 0; i<50 ; i++)
    {
      Rprintf("%d  ",  inedges[i].value);
      if ((i+1)%5 == 0)
	Rprintf("\n");
    }
#endif

  return g0;
}

/*  ***************** */
/*  void GraphDestroy */
/*  */
/*  ***************** */
void GraphDestroy0() {
  int i;

  for (i = 0; i < MAXEDGES; i++) {
    inedges0[i].value = outedges0[i].value = 0;  
    inedges0[i].parent = outedges0[i].parent = 0;  
    inedges0[i].left = outedges0[i].left = 0;  
    inedges0[i].right = outedges0[i].right = 0;  
  }

  for (i=0; i<=n_nodes; i++)
    outdegree0[i] = indegree0[i] = 0;

  free(outdegree0);
  free(indegree0);  
}
