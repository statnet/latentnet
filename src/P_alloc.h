#ifndef P_ALLOC_H
#define P_ALLOC_H

struct PMemNode_struct{
  void *data;
  struct PMemNode_struct *next;
};

typedef struct PMemNode_struct PMemNode;

void *P_alloc(size_t nmemb, size_t size);
void P_free_all();
void P_free_after(PMemNode *bookmark);
PMemNode *P_bookmark();
#endif
