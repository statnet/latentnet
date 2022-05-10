/*  File src/P_alloc.h in package latentnet, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
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
