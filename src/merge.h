#ifndef MERGE_H_INCLUDED
#define MERGE_H_INCLUDED

#include <stdio.h>   // for size_t (optional, kept for consistency)

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *try_merge_two_contigs(char *a, int len_a, char *b, int len_b, int max_subs, int min_overlap);
int contains_with_substitutions(const char *text, const char *pattern, int max_subs);
int merge_contigs(void);


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // MERGE_H_INCLUDED