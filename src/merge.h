#ifndef MERGE_H_INCLUDED
#define MERGE_H_INCLUDED

#include <stdio.h>   // for size_t (optional, kept for consistency)

typedef struct {
    int index_p;
    int index_s;
    int merged_length;
    int overlap_length;
    int substitutions;
} merge_info;

typedef struct {
    int * edges;
    int number_edges;
} edge_info;

typedef struct {
    int new_index_s;
    char* new_contig;
} contig_index_info;

typedef struct {
    char **final_contigs;
    int number_contigs;
} contig_array;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

merge_info try_merge_two_contigs(char *a, int len_a, char *b, int len_b, int min_overlap);
char* try_merge_two_reads(char *a, int len_a, char *b, int len_b, float max_subs, int min_overlap);
int contains_with_substitutions(const char *text, const char *pattern, int max_subs);
contig_array merge_contigs(int number_threads, int number_sequences);
char** filter_contained_contigs(char **contigs, int count, int *new_count);
contig_array select_merge(int number_threads, int number_sequences, int use_dependencies);


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // MERGE_H_INCLUDED