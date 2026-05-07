/**
 * recon.c - FCM-based read extension, contig assembly, and merging
 *
 * This program builds a Frequency Chaos Game Representation (FCGR) model
 * from a FASTQ file, extends reads that are not already covered by
 * existing contigs, then merges similar/overlapping contigs to reduce
 * redundancy.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include "hash.h"
#include "defs.h"
#include "merge.h"

// ----------------------------------------------------------------------
// Data structure for storing sequence positions
// ----------------------------------------------------------------------
typedef struct seq_info {
    short int used;
    long initial_position;
    long last_position;
} seq_info;

typedef struct read_pair {
    char *forward_read;
    char *reverse_read;
    short int valid; 
} read_pair;

// ----------------------------------------------------------------------
// Global options
// ----------------------------------------------------------------------
int number_of_threads = 1;
char *forward_file = NULL;
char *reverse_file = NULL;
char *output_path = "output.tsv";
int context_size;
int verbose = 0;
int help_menu = 0;

size_t capacity = 100000;  // initial capacity for sequence info array
seq_info *arr_forward;
seq_info *arr_reverse;
int last_sequence_id = 0;



HASH *hm;                      // hash map: kmer -> count

// ----------------------------------------------------------------------
// Contig assembly globals
// ----------------------------------------------------------------------
int max_substitutions = 3;                // substitution threshold (default)
int overlap_threshold = 30; // minimum overlap length for merging contigs

char **list_contigs_found = NULL;         // dynamic list of contigs
int num_contigs = 0;                      // number of contigs stored
int contig_capacity = 0;              // allocated size

// ----------------------------------------------------------------------
// Command line options
// ----------------------------------------------------------------------
static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"forward", required_argument, 0, 'f'},
    {"reverse", required_argument, 0, 'r'},
    {"output", required_argument, 0, 'o'},
    {"threads", required_argument, 0, 't'},
    {"context", required_argument, 0, 'c'},
    {"substitutions", required_argument, 0, 's'},
    {"overlap", required_argument, 0, 'e'},
    {"verbose", no_argument, 0, 'v'},
    {0, 0, 0, 0}
};

// ----------------------------------------------------------------------
// Helper: print usage
// ----------------------------------------------------------------------
void program_usage(const char *prog) {
    printf("\nUsage: %s -f <forward.fastq> -r <reverse.fastq> -c <context_size> [options]\n\n", prog);
    printf("Options:\n");
    printf("  -h, --help              Show this help\n");
    printf("  -f, --forward FILE      Forward FASTQ file (required)\n");
    printf("  -r, --reverse FILE      Reverse FASTQ file (required)\n");
    printf("  -o, --output FILE       Output TSV file (default: output.tsv)\n");
    printf("  -t, --threads N         Number of threads (default: 1)\n");
    printf("  -c, --context N         Context size (kmer length, >=2)\n");
    printf("  -s, --substitutions N   Max substitutions for merging (default: 3)\n");
    printf("  -e, --overlap N         Minimum overlap for merging contigs (default: 30)\n");
    printf("  -v, --verbose           Verbose output\n\n");
    help_menu = 1;
}

// ----------------------------------------------------------------------
// Parse command line arguments
// ----------------------------------------------------------------------
int option_parsing(int argc, char *argv[]) {
    int opt;
    if (argc <= 1) {
        program_usage(argv[0]);
        return 0;
    }
    while ((opt = getopt_long(argc, argv, "hf:r:o:t:c:vs:e:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'h': program_usage(argv[0]); return 0;
            case 'f': forward_file = optarg; break;
            case 'r': reverse_file = optarg; break;
            case 'o': output_path = optarg; break;
            case 't': number_of_threads = atoi(optarg); break;
            case 'c': context_size = atoi(optarg); break;
            case 'v': verbose = 1; break;
            case 's': max_substitutions = atoi(optarg); break;
            case 'e': overlap_threshold = atoi(optarg); break;
            default:  program_usage(argv[0]); return 1;
        }
    }
    if (!forward_file || !reverse_file || !context_size) {
        printf("Error: input files and context size are required.\n");
        return 1;
    }
    if (context_size < 2) {
        printf("Error: context size must be at least 2.\n");
        return 1;
    }
    if (access(output_path, F_OK) == 0) remove(output_path);
    return 0;
}

// ----------------------------------------------------------------------
// Hash table operations
// ----------------------------------------------------------------------
static void IncrementKmer(const char *kmer) { //Increment the count of a kmer in the hash table, or set to 1 if unseen
    int *val = GetValue(hm, kmer);
    if (val == NULL)
        SetValue(hm, kmer, 1);
    else
        SetValue(hm, kmer, *val + 1);
}

// ----------------------------------------------------------------------
// Dynamic array for sequence positions
// ----------------------------------------------------------------------
seq_info *expand_seq_info(seq_info *arr, size_t current_capacity, size_t extra) {
    size_t new_cap = current_capacity + extra;
    seq_info *new_arr = realloc(arr, new_cap * sizeof(seq_info));
    if (!new_arr) {
        fprintf(stderr, "realloc failed\n");
        free(arr);
        exit(1);
    }
    return new_arr;
}

void save_position(int id, long initial_position, long last_pos, seq_info *arr) {
    if (verbose)
        printf("Saving position for sequence %d: initial=%ld, last=%ld\n",
               id, initial_position, last_pos);
    if (id >= (int)capacity) {
        arr = expand_seq_info(arr, capacity, 200);
        capacity += 200;
    }
    arr[id].initial_position = initial_position;
    arr[id].last_position = last_pos;
    arr[id].used = 0;
}

// ----------------------------------------------------------------------
// Model training: read FASTQ, extract all kmers, update hash table
// ----------------------------------------------------------------------
void train_model(char *file_name, seq_info *arr) {
    FILE *file = fopen(file_name, "rb");
    if (!file) {
        perror("fopen");
        exit(1);
    }

    char buffer[BUFFER_SIZE];
    size_t bytes_read;
    long file_pos = 0;

    int part_sequence = 0;
    long initial_position = 0;
    long last_pos = 0;
    int id = 0;

    char *kmer_ring = malloc(context_size + 1);
    if (!kmer_ring) { perror("malloc"); exit(1); }
    int ring_pos = 0;
    int full_slots = 0;

    // Read file in chunks, process character by character
    while ((bytes_read = fread(buffer, 1, BUFFER_SIZE, file)) > 0) {
        for (size_t idx = 0; idx < bytes_read; ++idx) {
            
            char ch = buffer[idx];

            if (ch == '@') {
                if (id > 0)
                    save_position(id, initial_position, last_pos, arr); // save previous record's range
                part_sequence = 0;   // reset state
                ring_pos = 0;        // reset ring buffer
                full_slots = 0;      // reset k-mer window
                id++;                // increment record counter
            }
            else if (ch == '+') {
                part_sequence = 2; // quality header line, ignore
            }
            else if (ch == '\n') { 
                if (part_sequence == 0) { //end of header, save initial position of sequence
                    part_sequence = 1;
                    initial_position = file_pos + 1;
                }
            }
            else if (part_sequence == 1) { // in sequence part, process k-mers
                kmer_ring[ring_pos] = ch;
                ring_pos = (ring_pos + 1) % context_size;
                if (full_slots < context_size) full_slots++; //if the ring buffer isn't full yet, increment the count of filled slots
                if (full_slots == context_size) {
                    char kmer_str[context_size + 1];
                    int start = ring_pos;
                    for (int i = 0; i < context_size; ++i) { // build k-mer string from ring buffer
                        kmer_str[i] = kmer_ring[(start + i) % context_size];
                    }
                    kmer_str[context_size] = '\0';
                    IncrementKmer(kmer_str); // update hash table count
                }
                last_pos = file_pos;
            }
            file_pos++;
        }
    }
    if (id > 0) //save the last record's position
        save_position(id, initial_position, last_pos, arr);
    last_sequence_id = id;

    free(kmer_ring);
    fclose(file);
}

// ======================================================================
// Contig generation functions
// ======================================================================

char generate_next_symbol(const char *context, int forward) {
    char bases[] = {'A', 'T', 'C', 'G'};
    int counts[4] = {0};
    int max_count = 0;
    int best_indices[4];
    int num_best = 0;

    for (int i = 0; i < 4; ++i) {
        char key[context_size + 1];
        if (forward) {
            snprintf(key, sizeof(key), "%s%c", context, bases[i]);
        } else {
            snprintf(key, sizeof(key), "%c%s", bases[i], context);
        }
        int *val = GetValue(hm, key);
        counts[i] = (val == NULL) ? 0 : *val;
        if (counts[i] > max_count) {
            max_count = counts[i];
            num_best = 0;
            best_indices[num_best++] = i;          // new winner, reset list
        } else if (counts[i] == max_count && max_count > 0) {
            best_indices[num_best++] = i;          // tie, add to list
        }
    }

    if (max_count == 0) return '\0';                    // dead end
    int idx = best_indices[rand() % num_best];          // random tiebreak
    return bases[idx];
}

char * generate_sequence(const char *start_kmer, int forward) {
    char *context = malloc(context_size +1); // +1 for null terminator
    if (!context) { perror("malloc"); exit(1); }

    if (forward) {
        strncpy(context, start_kmer + 1, context_size - 1);
        context[context_size - 1] = '\0';
    } else {
        strncpy(context, start_kmer, context_size - 1);
        context[context_size - 1] = '\0';
    }

    int max_len = 4096; // initial max length for generated sequence
    char *sequence = malloc(max_len + 1);
    if (!sequence) { perror("malloc"); exit(1); }
    sequence[0] = '\0';
    int seq_len = 0;

    char next;
    while ((next = generate_next_symbol(context, forward)) != '\0') {
        if (seq_len + 1 >= max_len) { 
            return sequence; // reached max length, return what we have
        }

        sequence[seq_len++] = next;
        sequence[seq_len] = '\0';

        //printf("Generated next symbol '%c' with context '%s' -> sequence so far: %d\n", next, context, seq_len); 

        if (forward) { // shift context left and add next at the end
            memmove(context, context + 1, context_size - 2);
            context[context_size - 2] = next;
            context[context_size - 1] = '\0';
        } else { // shift context right and add next at the beginning
            memmove(context + 1, context, context_size - 2);
            context[0] = next;
            context[context_size - 1] = '\0';
        }
    }

    free(context);
    return sequence;
}


char *get_read_from_coordinates(long begin, long end, char *input_file) {
    FILE *f = fopen(input_file, "rb");
    if (!f) {
        perror("fopen in get_read_from_coordinates");
        return NULL;
    }
    fseek(f, begin, SEEK_SET);
    long len = end - begin + 1;
    char *read = malloc(len + 1);
    if (!read) { perror("malloc"); fclose(f); return NULL; }
    size_t n = fread(read, 1, len, f);
    read[n] = '\0';

    char *src = read, *dst = read;
    while (*src) {
        if (*src != '\n' && *src != '+') *dst++ = *src;
        src++;
    }
    *dst = '\0';
    fclose(f);
    return read;
}

char *extend_read(const char *read) {
    int read_len = strlen(read);
    if (read_len < context_size) { // Read is shorter than the context size, skipping extension
        if (verbose) fprintf(stderr, "Read too short, skipping extension\n");
        return NULL;
    }

    char first_kmer[context_size + 1];
    char last_kmer[context_size + 1];
    strncpy(first_kmer, read, context_size); // first k-mer from the read
    first_kmer[context_size] = '\0';
    strncpy(last_kmer, read + read_len - context_size, context_size); // last k-mer from the read
    last_kmer[context_size] = '\0';

    // Generate the following bases until we hit a dead end (no next symbol with count > 0)
    char *end_ext = generate_sequence(last_kmer, 1);
    char *begin_ext = generate_sequence(first_kmer, 0);

    int begin_len = strlen(begin_ext);
    for (int i = 0; i < begin_len / 2; ++i) { // reverse the backward extension to get the correct order
        char tmp = begin_ext[i];
        begin_ext[i] = begin_ext[begin_len - 1 - i];
        begin_ext[begin_len - 1 - i] = tmp;
    }

    /*if (verbose) {
        printf("\n--- Extending read (first 60 bases shown) ---\n");
        printf("Original read (first 60): %.60s\n", read);
        printf("Backward extension (%d bases): %s\n", begin_len, begin_ext);
        printf("Forward extension (%d bases): %s\n", (int)strlen(end_ext), end_ext);
        printf("---------------------------------------------------\n");
    }*/

    int total_len = begin_len + read_len + strlen(end_ext); // +1 for null terminator
    char *contig = malloc(total_len + 1);
    if (!contig) { perror("malloc"); free(begin_ext); free(end_ext); return NULL; }
    strcpy(contig, begin_ext);   // copy backward extension first
    strcat(contig, read);        // append original read
    strcat(contig, end_ext);     // append forward extension

    free(begin_ext);
    free(end_ext);

    return contig;
}

read_pair get_read_to_extend(int total_seqs, int extended_count, int seq_id) {

    if (total_seqs - extended_count <= 0) return (read_pair){0};

    int rand_id = rand() % (total_seqs - extended_count); // random ID in the range of remaining sequences
    int counter = 0;
    read_pair reads;

    while (counter < total_seqs) {
        if (arr_forward[counter].used == 0) {
            if (rand_id == 0) {
                arr_forward[counter].used = 1; // mark as used
                reads.forward_read = get_read_from_coordinates(arr_forward[counter].initial_position, arr_forward[counter].last_position, forward_file);
                reads.reverse_read = get_read_from_coordinates(arr_reverse[counter].initial_position, arr_reverse[counter].last_position, reverse_file);
                reads.valid = 1;
                return reads;

            }
            rand_id--;
        }
        counter++;
    }
    return (read_pair){0}; // should not reach here if counts are correct

}

/**
 * Compute the Hamming distance between two strings of given length.
 * Stops early if the distance exceeds max_allowed.
 * Returns the distance (or a value > max_allowed if exceeded).
 */
int hamming_distance_limit(const char *s1, const char *s2, int len, int max_allowed) {
    int dist = 0;
    for (int i = 0; i < len; i++) {
        if (s1[i] != s2[i]) {
            dist++;
            if (dist > max_allowed) {
                break;  // early exit, already too high
            }
        }
    }
    return dist;
}

void GenerateContigs(int number_sequences) {
    int extended_count = 0;
    while (extended_count < number_sequences) {
        read_pair r = get_read_to_extend(number_sequences, extended_count, 0);
        if (!r.valid) break;

        int found_forward = 0;
        for (int i = 0; i < num_contigs; ++i) {
            if (contains_with_substitutions(list_contigs_found[i], r.forward_read, max_substitutions)) {
                found_forward = 1;
                break;
            }
        }

        if (!found_forward) {

            // Extend the read in both directions
            char *extended_forward_read = extend_read(r.forward_read);
            char *extended_reverse_read = extend_read(r.reverse_read);

            // If either extension failed, we cannot use them.
            if (!extended_forward_read || !extended_reverse_read) {
                free(extended_forward_read);
                free(extended_reverse_read);
                free(r.forward_read);
                free(r.reverse_read);
                extended_count++;
                continue;
            }

            // Try to merge the two extended reads
            char *merged_contig = try_merge_two_contigs(
                extended_forward_read, strlen(extended_forward_read),
                extended_reverse_read, strlen(extended_reverse_read),
                 max_substitutions, overlap_threshold);
            

            // Decide which string to keep
            char *contig_to_add = merged_contig ? merged_contig : extended_forward_read;
            int substituted = 0;

            // Check if this contig replaces an existing one
            if (num_contigs > 0) {
                for (int i = 0; i < num_contigs && !substituted; ++i) {
                    if (contains_with_substitutions(contig_to_add, list_contigs_found[i], max_substitutions)) {
                        free(list_contigs_found[i]);          // free old contig
                        list_contigs_found[i] = contig_to_add; // replace with new
                        substituted = 1;
                        if (verbose) {
                            printf("Substituted contig %d with new contig of length %zu\n", i, strlen(contig_to_add));
                        }
                        break;
                    }
                }
            }

            if (!substituted || num_contigs == 0) { // If not substituted, add as new contig
                if (num_contigs >= contig_capacity) {
                    contig_capacity += 1000;
                    list_contigs_found = realloc(list_contigs_found, contig_capacity * sizeof(char*));
                    if (!list_contigs_found) { perror("realloc"); exit(1); }
                }
                list_contigs_found[num_contigs++] = contig_to_add;
                if (verbose) printf("Added contig of length %zu (total %d)\n", strlen(contig_to_add), num_contigs);
            }

            // Free the string that was NOT used
            if (merged_contig) {
                free(extended_forward_read);
                free(extended_reverse_read);
            } else {
                free(extended_reverse_read);
                // extended_forward_read is now owned by contig_to_add
            }
        } else if (verbose) {
            printf("Read already covered by existing contig\n");
        }

        free(r.forward_read);
        free(r.reverse_read);
        extended_count++;
    }
}

// ======================================================================
// Comparison function for qsort: sort by string length (ascending)
int compare_contig_len(const void *a, const void *b) {
    const char *sa = *(const char **)a;
    const char *sb = *(const char **)b;
    size_t len_a = strlen(sa);
    size_t len_b = strlen(sb);
    if (len_a < len_b) return 1;
    if (len_a > len_b) return -1;
    return 0;
}

// ======================================================================
// Main
// ======================================================================
int main(int argc, char *argv[]) {

    srand(time(NULL));

    // Menu and option parsing
    int rc = option_parsing(argc, argv);
    if (rc != 0) exit(rc);
    if (help_menu) exit(0);

    //Allocate sequence info array  - storing initial and last positions of each sequence in the FASTQ file
    arr_forward = malloc(capacity * sizeof(seq_info));
    if (!arr_forward) { perror("malloc"); exit(1); }
    arr_reverse = malloc(capacity * sizeof(seq_info));
    if (!arr_reverse) { perror("malloc"); exit(1); }


    if (access(forward_file, F_OK) != 0) {
        printf("Error: forward file does not exist.\n");
        exit(1);
    }
    if (access(reverse_file, F_OK) != 0) {
        printf("Error: reverse file does not exist.\n");
        exit(1);
    }

    printf("Running with %d thread(s)...\n", number_of_threads);
    printf("Training model with context size %d...\n", context_size);

    // Create hash table and train model on input FASTQ files
    hm = CreateHashTable();
    train_model(forward_file, arr_forward);
    train_model(reverse_file, arr_reverse);

    printf("\n--- Generating contigs ---\n");
    GenerateContigs(last_sequence_id);

    printf("\n--- Sort contigs by length ---\n");

    qsort(list_contigs_found, num_contigs, sizeof(char *), compare_contig_len);



    // Merge similar/overlapping contigs
    printf("\n--- Merging contigs ---\n");
    int merges_done;
    do {
        merges_done = merge_contigs();

        if (merges_done < 0) {
            fprintf(stderr, "Error during merging contigs.\n");
            break;
        } else {
            if (verbose) printf("Merged %d pairs, %d contigs remaining\n", merges_done, num_contigs);
        }
    } while (merges_done > 0);

    int total_contigs = num_contigs; // store total before output

    // Output final contigs
    printf("\n=== Final Contigs ===\n");
    for (int i = 0; i < num_contigs; ++i) {
        //printf(">contig_%d (length %zu)\n%s\n", i+1, strlen(list_contigs_found[i]), list_contigs_found[i]);
        free(list_contigs_found[i]);
    }
    free(list_contigs_found);

    printf("\nTotal contigs: %d     %d\n", total_contigs, num_contigs);

    RemoveHashTable(hm);
    free(arr_forward);
    free(arr_reverse);
    printf("Done.\n");
    return 0;
}