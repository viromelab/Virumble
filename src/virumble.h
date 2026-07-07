#ifndef VIRUMBLE_H
#define VIRUMBLE_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <pthread.h>

#include "hash.h"
#include "defs.h"
#include "merge.h"
#include "aux.h"

// ----------------------------------------------------------------------
// Data structures
// ----------------------------------------------------------------------
typedef struct {
    short int used;
    unsigned long long initial_position;
    unsigned long long last_position;
} seq_info;

typedef struct {
    char *forward_read;
    char *reverse_read;
    short int valid; 
} read_pair;

// Structure to hold all data a thread needs
typedef struct {
    int thread_id;
    int number_threads;
    int total_sequences;
    int file_index; // 0 for forward, 1 for reverse, 2 for additional, 3 for paired
} thread_data_t;

typedef struct {
    int thread_id;
    int number_threads;
    int * edges;
    int number_edges;
    merge_info *merge_possibilities;
    int number_merge_possibilities;
} thread_data_edges_t;

typedef struct {
    char *contig;
    int used;
} contig_info;

// ----------------------------------------------------------------------
// Global variables (extern, defined in virumble.c)
// ----------------------------------------------------------------------
extern int number_threads;
extern char *forward_file;
extern char *reverse_file;
extern char *additional_file;
extern char *output_path;
extern int context_size;
extern int verbose;
extern int help_menu;
extern size_t capacity;
extern seq_info *arr_forward;
extern seq_info *arr_reverse;
extern seq_info *arr_additional;
extern int last_sequence_id;
extern int seq_count;
extern int min_length;
extern pthread_mutex_t count_mutex;
extern pthread_mutex_t contig_mutex;
extern HASH *hm;
extern int max_substitutions;
extern int overlap_threshold;
extern contig_info *list_contigs;
extern int num_contigs;
extern int contig_capacity;

// ----------------------------------------------------------------------
// Function prototypes
// ----------------------------------------------------------------------
void program_usage(const char *prog);
int option_parsing(int argc, char *argv[]);
void IncrementKmer(const char *kmer);
seq_info *expand_seq_info(seq_info *arr, size_t current_capacity, size_t extra);
seq_info* save_position(int id, unsigned long long initial_position, unsigned long long last_pos, seq_info *arr, size_t *arr_capacity);
int train_model(char *file_name, seq_info **arr_ptr, size_t *cap_ptr);
char generate_next_symbol(const char *context, int forward);
char *generate_sequence(const char *start_kmer, int forward);
char *get_read_from_coordinates(unsigned long long begin, unsigned long long end, char *input_file);
char *extend_read(const char *read, int read_id, char *input_file);
read_pair get_read_to_extend(int total_seqs, int extended_count, int seq_id);
int hamming_distance_limit(const char *s1, const char *s2, int len, int max_allowed);
void update_count(void);
void add_contig_to_list(char *contig);
int compare_contig_len(const void *a, const void *b);
void output_contigs(char **contigs, int num_contigs, int min_length);

#endif // VIRUMBLE_H