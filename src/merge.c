#include "merge.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// ---------- External globals (defined in the main program) ----------
extern int verbose;
extern int max_substitutions;
extern int context_size;
extern char **list_contigs_found;
extern int num_contigs;
extern int contig_capacity;
extern int overlap_threshold;

// Heuristic: consider overlaps longer than this as unlikely to be useful.
// Set to 0 to disable limit (check all possible overlaps).
#define MAX_OVERLAP 50

// ---------- Helper: fast Hamming distance with early exit ----------
static int hamming_distance_limit(const char *s1, const char *s2, int len, int max_subs) {
    int dist = 0;
    for (int i = 0; i < len; ++i) {
        if (s1[i] != s2[i]) {
            if (++dist > max_subs) return dist;
        }
    }
    return dist;
}

// ---------- Try to merge two contigs (optimised) ----------
char *try_merge_two_contigs(char *a, int len_a, char *b, int len_b, int max_subs, int min_overlap) {
    // Overlap cannot exceed either length
    int max_possible = len_a < len_b ? len_a : len_b;
    if (min_overlap > max_possible) return NULL;

    // Limit overlap search to a reasonable maximum (speeds up long contigs)
    int overlap_start = max_possible;
    if (MAX_OVERLAP > 0 && overlap_start > MAX_OVERLAP)
        overlap_start = MAX_OVERLAP;
    // Search from large overlaps down to min_overlap (find best merge first)
    for (int overlap = overlap_start; overlap >= min_overlap; --overlap) {
        // Orientation 1: suffix of a vs prefix of b
        int dist = hamming_distance_limit(a + len_a - overlap, b, overlap, max_subs);
        if (dist <= max_subs) {
            int merged_len = len_a + len_b - overlap;
            char *merged = (char*)malloc(merged_len + 1);
            if (!merged) return NULL;
            memcpy(merged, a, len_a);
            memcpy(merged + len_a, b + overlap, len_b - overlap);
            merged[merged_len] = '\0';
            return merged;
        }

        // Orientation 2: suffix of b vs prefix of a
        dist = hamming_distance_limit(b + len_b - overlap, a, overlap, max_subs);
        if (dist <= max_subs) {
            int merged_len = len_b + len_a - overlap;
            char *merged = (char*)malloc(merged_len + 1);
            if (!merged) return NULL;
            memcpy(merged, b, len_b);
            memcpy(merged + len_b, a + overlap, len_a - overlap);
            merged[merged_len] = '\0';
            return merged;
        }
    }
    return NULL;
}

// ---------- Public functions ----------
int contains_with_substitutions(const char *text, const char *pattern, int max_subs) {
    int pattern_len = (int)strlen(pattern);
    int text_len = (int)strlen(text);
    if (pattern_len > text_len) return 0;

    // Fast path for exact matching
    if (max_subs == 0) {
        return strstr(text, pattern) != NULL;
    }

    // Sliding window with early exit
    for (int i = 0; i <= text_len - pattern_len; ++i) {
        int mismatches = 0;
        for (int j = 0; j < pattern_len; ++j) {
            if (text[i + j] != pattern[j]) {
                if (++mismatches > max_subs) break;
            }
        }
        if (mismatches <= max_subs) return 1;
    }
    return 0;
}

int merge_contigs(void) {
    if (num_contigs < 2) return -1;

    int merged_count = 0;

    // Precompute lengths once (will be updated when merging)
    int *lens = (int*)malloc(num_contigs * sizeof(int));
    if (!lens) { perror("malloc"); return -1; }
    for (int i = 0; i < num_contigs; ++i)
        lens[i] = (int)strlen(list_contigs_found[i]);

    int changed;
    do {
        changed = 0;

        // ----- Phase 1: remove contigs that are substrings of others (optional, disabled by default) -----
        // Uncomment if needed; now optimised with lengths
        
        int *to_remove = (int*)calloc(num_contigs, sizeof(int));
        if (!to_remove) { perror("calloc"); free(lens); return -1; }

        for (int i = 0; i < num_contigs; ++i) {
            if (to_remove[i]) continue;
            for (int j = 0; j < num_contigs; ++j) {
                if (i == j || to_remove[j]) continue;
                // Only check if lens[i] <= lens[j] (shorter can be inside longer)
                if (lens[i] <= lens[j] &&
                    contains_with_substitutions(list_contigs_found[j],
                                                list_contigs_found[i],
                                                max_substitutions)) {
                    to_remove[i] = 1;
                    changed = 1;
                    break;
                }
            }
        }

        // Compact and free removed strings
        int new_count = 0;
        for (int i = 0; i < num_contigs; ++i) {
            if (!to_remove[i]) {
                list_contigs_found[new_count] = list_contigs_found[i];
                lens[new_count] = lens[i];
                new_count++;
            } else {
                free(list_contigs_found[i]);
            }
        }
        num_contigs = new_count;
        free(to_remove);
        if (num_contigs < 2) break;
        

        // ----- Phase 2: merge overlapping pairs -----
        int merged_flag = 0;
        // Only consider unordered pairs (i < j) to avoid duplicate work
        for (int i = 0; i < num_contigs && !merged_flag; ++i) {
            for (int j = i + 1; j < num_contigs; ++j) {
                char *merged = try_merge_two_contigs(list_contigs_found[i], lens[i],
                                                     list_contigs_found[j], lens[j],
                                                     max_substitutions, overlap_threshold);
                if (merged) {
                    if (verbose) {
                        printf("Merging contig %d (len %d) and %d (len %d) -> new length %zu\n",
                               i, lens[i], j, lens[j], strlen(merged));
                    }
                    merged_count++;
                    free(list_contigs_found[i]);
                    list_contigs_found[i] = merged;
                    lens[i] = (int)strlen(merged);
                    free(list_contigs_found[j]);

                    // Remove contig j by shifting later entries left
                    for (int k = j; k < num_contigs - 1; ++k) {
                        list_contigs_found[k] = list_contigs_found[k + 1];
                        lens[k] = lens[k + 1];
                    }
                    num_contigs--;
                    merged_flag = 1;
                    changed = 1;
                    break;   // restart from the beginning after a successful merge
                }
            }
        }
    } while (changed);

    free(lens);

    // Shrink capacity if needed
    if (num_contigs < contig_capacity) {
        list_contigs_found = (char**)realloc(list_contigs_found,
                                             num_contigs * sizeof(char*));
        if (!list_contigs_found && num_contigs > 0) {
            perror("realloc");
            exit(1);
        }
        contig_capacity = num_contigs;
    }
    return merged_count;
}