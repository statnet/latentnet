#ifndef P_ALLOC_H
#define P_ALLOC_H

// If DEBUG is not set, just use R_alloc.
#ifndef DEBUG
#define P_alloc(nmemb, size) R_alloc(nmemb, size)
#define P_free_all()
#define P_free_after(bookmark)
#define P_bookmark()
#endif

struct PMemNode_struct{
  void *data;
  struct PMemNode_struct *next;
};

typedef struct PMemNode_struct PMemNode;

#ifdef DEBUG
void *P_alloc(size_t nmemb, size_t size);
void P_free_all();
void P_free_after(PMemNode *bookmark);
PMemNode *P_bookmark();
#endif

#endif
