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
  printf("%p",PMemNodes);
  if(PMemNodes) printf("(%p)",PMemNodes->data);
}

void *P_alloc(size_t nmemb, size_t size){
  printf("P_alloc: ");
  P_print_alloc();
  printf(" -%uB-> ", nmemb*size);

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
  printf("\n");
  return memnode->data;
}

/*R_INLINE*/ void P_free_after(PMemNode *bookmark){
  while(PMemNodes!=bookmark && PMemNodes){

    printf("P_free: ");
    P_print_alloc();

    free(PMemNodes->data);
    PMemNode *temp=PMemNodes;
    PMemNodes=PMemNodes->next;
    free(temp);

    printf(" --> ");
    P_print_alloc();
    printf("\n");
  }
}

void P_free_all(){
  P_free_after(NULL);
}

PMemNode *P_bookmark(){
  return PMemNodes;
}

#endif /*DEBUG*/
