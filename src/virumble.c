/**
 * virumble.c - FCM-based read extension, contig assembly, and merging
 *
 * This program builds a Frequency Chaos Game Representation (FCGR) model
 * from a FASTQ file, extends reads that are not already covered by
 * existing contigs, then merges similar/overlapping contigs to reduce
 * redundancy.
 */

#include "virumble.h"

// ----------------------------------------------------------------------
// Global options (definitions)
// ----------------------------------------------------------------------
int number_of_threads = 1;

char *forward_file = NULL;
char *reverse_file = NULL;
char *additional_file = NULL;

char *output_path = "output.fa";

int context_size;
int verbose = 0;
int help_menu = 0;

size_t capacity = 10000;  // initial capacity for sequence info array
seq_info *arr_forward;
seq_info *arr_reverse;
seq_info *arr_additional; // only used if additional file is provided
int last_sequence_id = 0;
int seq_count = 0; // count of sequences processed for extension
int min_length = 0; // minimum length of contigs to output


int num_reads = 0;
int length_all_sequences = 0;

int avg_read_length = 0;


pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t contig_mutex = PTHREAD_MUTEX_INITIALIZER;

HASH *hm;                      // hash map: kmer -> count

// ----------------------------------------------------------------------
// Contig assembly globals
// ----------------------------------------------------------------------
int max_overlap = -1;
int min_overlap = -1;              // minimum overlap length for merging contigs
int max_substitutions = -1;                // substitution threshold (default)

float f_max_substitutions = 0.05;
float f_max_overlap = 1.0;
float f_min_overlap = 0.1;

contig_info *list_contigs = NULL;        // dynamic list of contigs
int num_contigs = 0;                     // number of contigs stored
int contig_capacity = 0;                 // allocated size

// ----------------------------------------------------------------------
// Command line options
// ----------------------------------------------------------------------
static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"forward", required_argument, 0, 'f'},
    {"reverse", required_argument, 0, 'r'},
    {"additional", required_argument, 0, 'a'},
    {"output", required_argument, 0, 'o'},
    {"context", required_argument, 0, 'c'},
    {"substitutions", required_argument, 0, 's'},
    {"max_overlap", required_argument, 0, 'M'},
    {"min_overlap", required_argument, 0, 'm'},
    {"threads", required_argument, 0, 't'},
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
    printf("  -a, --additional FILE   Additional FASTQ file (optional)\n");
    printf("  -o, --output FILE       Output FASTA file (default: output.fa)\n");
    printf("  -c, --context N         Context size (kmer length, >=2)\n");
    printf("  -s, --substitutions N   Maximum number of substitutions for merging (can be set as the number of bases or as a percentage of the average read length)\n");
    printf("  -M, --max_overlap N     Maximum overlap for merging contigs (can be set as the number of bases or as a percentage of the average read length)\n");
    printf("  -m, --min_overlap N     Minimum overlap for merging contigs (can be set as the number of bases or as a percentage of the average read length)\n");
    printf("  -l, --length N          Minimum length for output contigs (default: 0)\n");
    printf("  -t, --threads N         Number of threads (default: 1)\n");
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
    while ((opt = getopt_long(argc, argv, "hf:r:a:o:c:s:M:m:l:t:v", long_options, NULL)) != -1) {
        switch (opt) {
            case 'h': program_usage(argv[0]); return 0;
            case 'f': forward_file = optarg; break;
            case 'r': reverse_file = optarg; break;
            case 'a': additional_file = optarg; break;
            case 'o': output_path = optarg; break;
            case 'c': context_size = atoi(optarg); break;

            case 's': {
                double val;
                NumType type = parse_number(optarg, &val);
                if (type == TYPE_INT) {
                    max_substitutions = val;
                } else if (type == TYPE_FLOAT) {
                    f_max_substitutions = val;
                } else {
                    fprintf(stderr, "Error: -s requires a valid number (integer or float)\n");
                    program_usage(argv[0]);
                    return 1;
                }
                break;
            }
            case 'M': {
                double val;
                NumType type = parse_number(optarg, &val);
                if (type == TYPE_INT) {
                    max_overlap = val; 
                } else if (type == TYPE_FLOAT) {
                    f_max_overlap = val;
                } else {
                    fprintf(stderr, "Error: -M requires a valid number (integer or float)\n");
                    program_usage(argv[0]);
                    return 1;
                }
                break;
            }
            case 'm': {
                double val;
                NumType type = parse_number(optarg, &val);
                if (type == TYPE_INT) {
                    min_overlap = val;
                } else if (type == TYPE_FLOAT) {
                    f_min_overlap = val;
                } else {
                    fprintf(stderr, "Error: -m requires a valid number (integer or float)\n");
                    program_usage(argv[0]);
                    return 1;
                }
                break;
            }

            case 'l': min_length = atoi(optarg); break;
            case 't': number_of_threads = atoi(optarg); break;
            case 'v': verbose = 1; break;

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

    

    
    // TODO - chaneg it so max_overlap and min_overlap are checked !!!



    if (access(output_path, F_OK) == 0) remove(output_path);

    return 0;
}

// ----------------------------------------------------------------------
// Hash table operations
// ----------------------------------------------------------------------
void IncrementKmer(const char *kmer) {
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
int train_model(char *file_name, seq_info *arr) {
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
    int id = -1;

    char *kmer_ring = malloc(context_size + 1);
    if (!kmer_ring) { perror("malloc"); exit(1); }
    int ring_pos = 0;
    int full_slots = 0;
    int length_sequences_curr_file = 0;

    while ((bytes_read = fread(buffer, 1, BUFFER_SIZE, file)) > 0) {
        for (size_t idx = 0; idx < bytes_read; ++idx) {
            char ch = buffer[idx];

            if (ch == '@') {
                if (id >= 0)
                    save_position(id, initial_position, last_pos, arr);
                //Reset values
                part_sequence = 0;
                ring_pos = 0;
                full_slots = 0;
                id++;
            }
            else if (ch == '+') {
                part_sequence = 2;
            }
            else if (ch == '\n') { 
                if (part_sequence == 0) {
                    part_sequence = 1;
                    initial_position = file_pos + 1;
                }
            }
            else if (part_sequence == 1) {
                kmer_ring[ring_pos] = ch;
                ring_pos = (ring_pos + 1) % context_size;
                if (full_slots < context_size) full_slots++;
                if (full_slots == context_size) {
                    char kmer_str[context_size + 1];
                    int start = ring_pos;
                    for (int i = 0; i < context_size; ++i) {
                        kmer_str[i] = kmer_ring[(start + i) % context_size];
                    }
                    kmer_str[context_size] = '\0';
                    IncrementKmer(kmer_str);
                }
                length_sequences_curr_file++;
                last_pos = file_pos;
            }
            file_pos++;
        }
    }
    if (id >= 0)
        save_position(id, initial_position, last_pos, arr);

    free(kmer_ring);
    fclose(file);

    num_reads += id + 1;
    length_all_sequences += length_sequences_curr_file;

    return id;
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
            best_indices[num_best++] = i;
        } else if (counts[i] == max_count && max_count > 0) {
            best_indices[num_best++] = i;
        }
    }

    if (max_count == 0) return '\0';
    int idx = best_indices[rand() % num_best];
    return bases[idx];
}

char *generate_sequence(const char *start_kmer, int forward) {
    char *context = malloc(context_size + 1);
    if (!context) { perror("malloc"); exit(1); }

    if (forward) {
        strncpy(context, start_kmer + 1, context_size - 1);
        context[context_size - 1] = '\0';
    } else {
        strncpy(context, start_kmer, context_size - 1);
        context[context_size - 1] = '\0';
    }

    int max_len = 4096;
    char *sequence = malloc(max_len + 1);
    if (!sequence) { perror("malloc"); exit(1); }
    sequence[0] = '\0';
    int seq_len = 0;

    char next;
    while ((next = generate_next_symbol(context, forward)) != '\0') {
        if (seq_len + 1 >= max_len) {
            return sequence;
        }
        sequence[seq_len++] = next;
        sequence[seq_len] = '\0';

        if (forward) {
            memmove(context, context + 1, context_size - 2);
            context[context_size - 2] = next;
            context[context_size - 1] = '\0';
        } else {
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

char *extend_read(const char *read, int read_id, char *input_file) {
    int read_len = strlen(read);
    if (read_len < context_size) {
        if (verbose) fprintf(stderr, "Read too short (%d bp), skipping extension (ID: %d ; File: %s) ---%s---\n", read_len, read_id, input_file, read);
        return NULL;
    }

    char first_kmer[context_size + 1];
    char last_kmer[context_size + 1];
    strncpy(first_kmer, read, context_size);
    first_kmer[context_size] = '\0';
    strncpy(last_kmer, read + read_len - context_size, context_size);
    last_kmer[context_size] = '\0';

    char *end_ext = generate_sequence(last_kmer, 1);
    char *begin_ext = generate_sequence(first_kmer, 0);

    int begin_len = strlen(begin_ext);
    for (int i = 0; i < begin_len / 2; ++i) {
        char tmp = begin_ext[i];
        begin_ext[i] = begin_ext[begin_len - 1 - i];
        begin_ext[begin_len - 1 - i] = tmp;
    }

    int total_len = begin_len + read_len + strlen(end_ext);
    char *contig = malloc(total_len + 1);
    if (!contig) { perror("malloc"); free(begin_ext); free(end_ext); return NULL; }
    strcpy(contig, begin_ext);
    strcat(contig, read);
    strcat(contig, end_ext);

    free(begin_ext);
    free(end_ext);

    return contig;
}


void update_count() {
    pthread_mutex_lock(&count_mutex);
    seq_count++;
    pthread_mutex_unlock(&count_mutex);
}

void add_contig_to_list(char *contig) {
    pthread_mutex_lock(&contig_mutex);
    if (num_contigs >= contig_capacity) {
        int new_capacity = contig_capacity == 0 ? 100 : contig_capacity * 2;
        contig_info *new_list = realloc(list_contigs, new_capacity * sizeof(contig_info));
        if (!new_list) {
            perror("realloc");
            exit(1);
        }
        list_contigs = new_list;
        contig_capacity = new_capacity;
    }
    list_contigs[num_contigs].contig = contig;
    list_contigs[num_contigs].used = 0;
    num_contigs++;
    pthread_mutex_unlock(&contig_mutex);
}

void *extend_all_reads(void *arg) {
    thread_data_t *data = (thread_data_t*) arg;
    int thread_id = data->thread_id;
    int total_sequences = data->total_sequences;
    int file_index = data->file_index;

    if (verbose) printf("Thread %d: Starting extension of single reads, file_index %d...\n", thread_id, file_index);

    seq_info *arr1;
    char *file1;

    seq_info *arr2 = NULL;
    char *file2 = NULL;

    switch (file_index) {
        case 0:
            arr1 = arr_forward;
            file1 = forward_file;
            break;
        case 1:
            arr1 = arr_reverse;
            file1 = reverse_file;
            break;
        case 2:
            arr1 = arr_additional;
            file1 = additional_file;
            break;
        case 3:
            arr1 = arr_forward;
            file1 = forward_file;
            arr2 = arr_reverse;
            file2 = reverse_file;
            break;
        default:
            printf("Invalid file index.\n");
            exit(1);
    }

    for (int i = 0; i < total_sequences; i++) {
        if (i % number_of_threads != thread_id) continue;

        update_count();

        char* read1 = get_read_from_coordinates(arr1[i].initial_position, arr1[i].last_position, file1);

        if (file_index == 3) {
            char* read2 = get_read_from_coordinates(arr2[i].initial_position, arr2[i].last_position, file2);

            // For this merge, no substitutions are acceptable and the overlap should be at least 8
            char *merged = try_merge_two_reads(read1, strlen(read1), read2, strlen(read2), 0.0, 8);

            if (merged) {
                printf("Thread %d: Merging paired-end reads of length %zu and %zu, id %d -> merged length %zu\n",
                       thread_id, strlen(read1), strlen(read2), i+1, strlen(merged));
                char *extended = extend_read(merged, i+1, "both files");

                if (extended == NULL) {
                    if (verbose) fprintf(stderr, "Thread %d: No extension possible for merged read of id %d, file %d\n", thread_id, i+1, file_index);
                }
                if (extended) add_contig_to_list(extended);
                free(merged);
            } else {
                char *extended_read1 = extend_read(read1, i+1, file1);
                if (extended_read1 == NULL) {
                    if (verbose) fprintf(stderr, "Thread %d: No extension possible for merged read of id %d, file %d\n", thread_id, i+1, file_index);
                }
                char *extended_read2 = extend_read(read2, i+1, file2);
                if (extended_read2 == NULL) {
                    if (verbose) fprintf(stderr, "Thread %d: No extension possible for merged read of id %d, file %d\n", thread_id, i+1, file_index);
                }
                char *merged_extended = try_merge_two_reads(read1, strlen(read1), read2, strlen(read2), f_max_substitutions, min_overlap);
                if (merged_extended) {
                    add_contig_to_list(merged_extended);
                    free(merged_extended);
                } else {
                    if (extended_read1) add_contig_to_list(extended_read1);
                    if (extended_read2) add_contig_to_list(extended_read2);
                }
            }
            free(read1);
            free(read2);
        } else {
            char *extended = extend_read(read1, i+1, "smt else");
            if (extended == NULL) {
                if (verbose) fprintf(stderr, "Thread %d: No extension possible for read of id %d, file %d\n", thread_id, i+1, file_index);
            }
            if (extended) add_contig_to_list(extended);
            free(read1);
        }
    }
    return NULL;
}

void GenerateContigs(int file_index, int number_sequences) {
    pthread_t threads[number_of_threads];
    thread_data_t thread_data[number_of_threads];

    for (int i = 0; i < number_of_threads; i++) {
        if (verbose) printf("Main: Starting thread %d\n", i);
        thread_data[i].thread_id = i;
        thread_data[i].total_sequences = number_sequences;
        thread_data[i].file_index = file_index;

        if (pthread_create(&threads[i], NULL, extend_all_reads, &thread_data[i]) != 0) {
            perror("Failed to create thread");
            exit(1);
        }
    }
    for (int i = 0; i < number_of_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    printf("%d sequences processed - %d contigs generated\n", seq_count, num_contigs);
    seq_count = 0;
}

// Comparison function for qsort: sort by length descending
int compare_contig_len(const void *a, const void *b) {
    const contig_info *ca = (const contig_info *)a;
    const contig_info *cb = (const contig_info *)b;
    size_t len_a = strlen(ca->contig);
    size_t len_b = strlen(cb->contig);
    if (len_a < len_b) return 1;
    if (len_a > len_b) return -1;
    return 0;
}

void output_contigs(char **contigs, int num_contigs, int min_length) {
    FILE *out = fopen(output_path, "a");
    if (!out) {
        perror("fopen output file");
        exit(1);
    }
    for (int i = 0; i < num_contigs; ++i) {
        size_t len = strlen(contigs[i]);
        if (len < (size_t)min_length) continue;
        //printf(">contig_%d (length %zu)\n%.60s\n", i+1, len, contigs[i]);
        fprintf(out, ">contig_%d (length %zu)\n%s\n", i+1, len, contigs[i]);
    }
    fclose(out);
}

void fix_overlap_values () {

    int avg_read_length = length_all_sequences / num_reads;

    if (max_substitutions == -1) {
        max_substitutions = f_max_substitutions * avg_read_length;
    } else {
        f_max_substitutions = (float) max_substitutions / avg_read_length;
        printf("%f    --- val max subs\n\n", f_max_substitutions);
    }

    if (max_overlap == -1) {
        max_overlap = f_max_overlap * avg_read_length;
    }

    if (min_overlap == -1) {
        min_overlap = f_min_overlap * avg_read_length;
    }

    if (min_overlap > max_overlap) {
        printf ("The minimum overlap is greater than the maximum overlap (Min overlap: %d  Max overlap:%d)\nExiting...\n", min_overlap, max_overlap);
        exit(1);
    }

}

// ======================================================================
// Main
// ======================================================================
int main(int argc, char *argv[]) {
    srand(time(NULL));

    int rc = option_parsing(argc, argv);
    if (rc != 0) exit(rc);
    if (help_menu) exit(0);

    if (access(forward_file, F_OK) != 0 || access(reverse_file, F_OK) != 0) {
        printf("Error: invalid input. Exiting...\n");
        exit(1);
    }
    if (additional_file && access(additional_file, F_OK) != 0) {
        printf("Error: invalid additional file. Exiting...\n");
        exit(1);
    }

    printf("Input files: %s, %s%s%s\n", forward_file, reverse_file,
           additional_file ? ", " : "", additional_file ? additional_file : "");
    printf("Training model with context size %d...\n", context_size);
    printf("Running with %d thread(s)...\n", number_of_threads);
    printf("Output path: %s\n", output_path);
    printf("Minimum length for output contigs: %d\n", min_length);
    printf("==============================\n");

    arr_forward = malloc(capacity * sizeof(seq_info));
    if (!arr_forward) { perror("malloc"); exit(1); }
    arr_reverse = malloc(capacity * sizeof(seq_info));
    if (!arr_reverse) { perror("malloc"); exit(1); }
    if (additional_file) {
        arr_additional = malloc(capacity * sizeof(seq_info));
        if (!arr_additional) { perror("malloc"); exit(1); }
    }

    printf("\n---Training model---\n");
    hm = CreateHashTable();

    int max_id_forward = train_model(forward_file, arr_forward);
    int max_id_reverse = train_model(reverse_file, arr_reverse);
    int max_id_additional = 0;
    if (additional_file) {
        max_id_additional = train_model(additional_file, arr_additional);
    }

    fix_overlap_values();

    printf("Max substitutions for merging: %d\n", max_substitutions);
    printf("Minimum overlap for merging contigs: %d\n", min_overlap);
    printf("Max overlap for merging: %d\n", max_overlap);

    printf("\n--- Generating contigs ---\n");
    char error_msg[256];
    int sync_result = check_paired_end_files(forward_file, reverse_file, error_msg, sizeof(error_msg));
    if (sync_result == 0) {
        printf("Paired-end files are synchronized.\n");
        GenerateContigs(3, max_id_forward);
    } else {
        fprintf(stderr, "Warning: Paired-end files are not synchronized. Reason: %s\nProceeding with reconstruction treating them as independent files.\n", error_msg);
        GenerateContigs(0, max_id_forward);
        GenerateContigs(1, max_id_reverse);
    }

    if (additional_file) {
        printf("\n--- Processing additional file ---\n");
        GenerateContigs(2, max_id_additional);
    }

    // Sort contigs by length descending (optional, uncomment if needed)
    // qsort(list_contigs, num_contigs, sizeof(contig_info), compare_contig_len);

    printf("\n--- Merging contigs ---\n");
    int total_contigs = num_contigs;
    contig_array merges_done = merge_contigs(number_of_threads, total_contigs);


    printf("\n=== Final Contigs ===\n");
    qsort(merges_done.final_contigs, merges_done.number_contigs, sizeof(char*), compare_contig_len);

    printf("Total contigs after merging: %d\n", merges_done.number_contigs);
    output_contigs(merges_done.final_contigs, merges_done.number_contigs, min_length);

    // Clean up

    free(list_contigs);
    RemoveHashTable(hm);
    free(arr_forward);
    free(arr_reverse);
    if (additional_file) free(arr_additional);
    printf("Done.\n");
    return 0;
}