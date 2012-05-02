/**********************************************************************/
/* Routines for memory management similar in functionality to R_alloc */
/* Used because it's easier to debug memory problems this way.        */
/* Turned off (using R_alloc instead) unless DEBUG is set.            */
/**********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "P_alloc.h"
#include <R.h>

PMemNode *PMemNodes=NULL;
#ifdef DEBUG
/*R_INLINE*/ void P_print_alloc(){
  Rprintf("%p",(void*) PMemNodes);
  if(PMemNodes) Rprintf("(%p)",PMemNodes->data);
}

void *P_alloc(size_t nmemb, size_t size){
  Rprintf("P_alloc: ");
  P_print_alloc();
  Rprintf(" --%uB-> ", (unsigned int) nmemb*size);

  PMemNode *memnode=(PMemNode *)calloc(1,sizeof(PMemNode));
  if(!memnode) return NULL;

  memnode->data=calloc(nmemb,size);
  if(!memnode->data){
    free(memnode);
    return NULL;
  }
  memnode->next=PMemNodes;
  PMemNodes=memnode;

  P_print_alloc();
  Rprintf("\n");
  return memnode->data;
}

/*R_INLINE*/ void P_free_after(PMemNode *bookmark){
  while(PMemNodes!=bookmark && PMemNodes){

    Rprintf("P_free: ");
    P_print_alloc();

    free(PMemNodes->data);
    PMemNode *temp=PMemNodes;
    PMemNodes=PMemNodes->next;
    free(temp);

    Rprintf(" --> ");
    P_print_alloc();
    Rprintf("\n");
  }
}

void P_free_all(){
  P_free_after(NULL);
}

PMemNode *P_bookmark(){
  return PMemNodes;
}

#endif /*DEBUG*/
