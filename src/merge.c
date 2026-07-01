#include "merge.h"
#include "virumble.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>
#include <float.h>
#include <limits.h>

// ---------- External globals (defined in the main program) ----------
extern int verbose;               // if >0, print detailed progress messages
extern int max_substitutions;     // maximum allowed mismatches when merging two contigs
extern float f_max_substitutions;
extern int avg_read_length;
extern contig_info *list_contigs; // array of input contigs (each with a string and metadata)
extern int min_overlap;     // minimum overlap length required to consider a merge
extern int number_threads;        // number of parallel threads to use
extern int max_overlap;           // maximum overlap length to consider (0 = no limit)
extern float f_max_substitutions_dedup;
extern int min_length;            // minimum length of contigs to output

pthread_mutex_t merge_mutex = PTHREAD_MUTEX_INITIALIZER; // protects the global merge list

merge_info *merge_possibilities = NULL; // dynamic array of possible merges found
int index_merge_possibilities   = 0;    // current number of stored merges
int merge_capacity              = 0;    // allocated size of merge_possibilities

// ---------- Helper: fast Hamming distance with early exit ----------
// Returns the number of mismatches between s1 and s2, but stops counting if it exceeds max_subs.
int hamming_distance_limit(const char *s1, const char *s2, int len, int max_subs) {
    int dist = 0;
    for (int i = 0; i < len; ++i) {
        if (s1[i] != s2[i]) {
            if (++dist > max_subs) return dist; // early exit: already too many mismatches
        }
    }
    return dist;
}

// ---------- Pairwise overlap check (no merge string produced) ----------
// Tests whether contigs a (len_a) and b (len_b) can be merged with at most max_subs mismatches
// over an overlap of at least min_overlap. Returns a merge_info structure describing the best merge.
// merge_info.index_p indicates which orientation: 1 = a as prefix, b as suffix; 2 = b as prefix, a as suffix.
// The actual indices ( index_p, index_s) are filled later by the caller.
merge_info try_merge_two_contigs(char *a, int len_a, char *b, int len_b, int min_overlap) {
    int max_possible = len_a < len_b ? len_a : len_b;
    if (min_overlap > max_possible) return (merge_info){0, 0, 0, 0, 0};

    int overlap_start = max_possible;
    if (max_overlap > 0 && overlap_start > max_overlap)  // enforce user‑supplied maximum overlap
        overlap_start = max_overlap;

    // Try decreasing overlap lengths (largest first gives the longest merged sequence)
    for (int overlap = overlap_start; overlap >= min_overlap; --overlap) {

        int max_subs = f_max_substitutions * avg_read_length;
        
        // Case 1: suffix of a overlaps prefix of b
        int dist = hamming_distance_limit(a + len_a - overlap, b, overlap, max_subs);
        if (dist <= max_subs) {
            int merged_len = len_a + len_b - overlap;
            return (merge_info){1, 0, merged_len, overlap, dist};
        }
        // Case 2: suffix of b overlaps prefix of a
        dist = hamming_distance_limit(b + len_b - overlap, a, overlap, max_subs);
        if (dist <= max_subs) {
            int merged_len = len_b + len_a - overlap;
            return (merge_info){2, 0, merged_len, overlap, dist};
        }
    }
    return (merge_info){0, 0, 0, 0, 0};
}

// ---------- Public utility ----------
// Checks whether 'pattern' appears as a substring of 'text' allowing up to 'max_subs' mismatches.
// If max_subs == 0, falls back to standard strstr.
int contains_with_substitutions(const char *text, const char *pattern, int max_subs) {
    int pattern_len = (int)strlen(pattern);
    int text_len    = (int)strlen(text);
    if (pattern_len > text_len) return 0;
    if (max_subs == 0) return strstr(text, pattern) != NULL;

    // Slide a window of length pattern_len over text
    for (int i = 0; i <= text_len - pattern_len; ++i) {
        int mismatches = 0;
        for (int j = 0; j < pattern_len; ++j) {
            if (text[i + j] != pattern[j])
                if (++mismatches > max_subs) break;
        }
        if (mismatches <= max_subs) return 1;
    }
    return 0;
}

// ---------- Thread worker: collect pairwise overlaps ----------
// Thread‑safe addition of a merge_info to the global list.
static void add_merge_to_list(merge_info info) {
    if (verbose)
        printf("Merge possible: contig %d and contig %d -> %d bp, overlap %d, subs %d\n",
               info.index_p, info.index_s, info.merged_length,
               info.overlap_length, info.substitutions);

    pthread_mutex_lock(&merge_mutex);
    merge_possibilities[index_merge_possibilities++] = info;
    pthread_mutex_unlock(&merge_mutex);
}

// Thread worker function: each thread examines a subset of contig pairs (i<j) where i % number_threads == thread_id.
void *get_merge_possibilities(void *arg) {
    thread_data_t *data        = (thread_data_t *)arg;
    int thread_id              = data->thread_id;
    int total_sequences        = data->total_sequences;
    int number_threads         = data->number_threads;
    int64_t pairs_checked = 0;
    int     merges_found  = 0;

    // Iterate over all unordered pairs (i,j) with i<j, but only those where i % number_threads == thread_id
    for (int i = 0; i < total_sequences; i++) {
        if (i % number_threads != thread_id) continue;

        for (int j = i + 1; j < total_sequences; j++) {
            pairs_checked++;
            char *read1 = list_contigs[i].contig;
            char *read2 = list_contigs[j].contig;
            if (!read1 || !read2) continue;

            merge_info info = try_merge_two_contigs(read1, (int)strlen(read1),
                                                    read2, (int)strlen(read2),
                                                    min_overlap);
            if (info.merged_length > 0) {
                // Store the original contig indices (the smaller index in index_p, larger in index_s)
                if (info.index_p == 1) {
                    info.index_p = i;
                    info.index_s = j;
                } else {
                    info.index_p = j;
                    info.index_s = i;
                }
                merges_found++;
                add_merge_to_list(info);
            }
        }
    }
    if (verbose)
        printf("Thread %d: checked %lld pairs, found %d merges\n",
               thread_id, (long long)pairs_checked, merges_found);
    return NULL;
}

// =============================================================================
// MSA ENGINE (Multiple Sequence Alignment)
// =============================================================================
//
// Progressive MSA pipeline:
//   1. Pairwise distance matrix derived from overlap substitution rates
//   2. UPGMA guide tree
//   3. Profile-vs-profile Needleman-Wunsch alignment
//   4. Majority-vote consensus from final alignment columns
//
// Alphabet encoding: A=0, C=1, G=2, T=3, N=4, gap=5
// =============================================================================

#define MSA_ALPHA 6          // number of symbols: A,C,G,T,N,-
#define GAP_CHAR  '-'        // gap character

// Convert a base character to an index (0‑4 for A/C/G/T/N, 5 for gap)
static int base_index(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        case '-':           return 5;
        default:            return 4;   // treat everything else as N
    }
}

// Reverse mapping: index -> character
static char index_base(int i) {
    return "ACGTN-"[i];
}

// ---------- Profile structure ----------
// A profile is a position‑specific frequency matrix over the alphabet.
// It represents a set of aligned sequences.
typedef struct {
    float **freq;   // [ncols][MSA_ALPHA] – frequency of each base at each column
    int    ncols;   // number of columns in the profile
    int    nseqs;   // number of sequences that went into this profile
} Profile;

// Allocate a Profile with ncols columns (frequencies initially zero)
static Profile *profile_alloc(int ncols) {
    Profile *p = (Profile *)malloc(sizeof(Profile));
    if (!p) return NULL;
    p->ncols = ncols;
    p->nseqs = 0;
    p->freq  = (float **)malloc(ncols * sizeof(float *));
    if (!p->freq) { free(p); return NULL; }
    for (int i = 0; i < ncols; i++) {
        p->freq[i] = (float *)calloc(MSA_ALPHA, sizeof(float));
        if (!p->freq[i]) {
            for (int k = 0; k < i; k++) free(p->freq[k]);
            free(p->freq); free(p); return NULL;
        }
    }
    return p;
}

// Free all memory associated with a Profile
static void profile_free(Profile *p) {
    if (!p) return;
    for (int i = 0; i < p->ncols; i++) free(p->freq[i]);
    free(p->freq);
    free(p);
}

// Build a Profile from an array of aligned sequences (all of length ncols)
static Profile *profile_from_seqs(char **seqs, int nseqs, int ncols) {
    Profile *p = profile_alloc(ncols);
    if (!p) return NULL;
    p->nseqs = nseqs;
    // Count occurrences
    for (int s = 0; s < nseqs; s++)
        for (int c = 0; c < ncols; c++)
            p->freq[c][base_index(seqs[s][c])] += 1.0f;
    // Convert to frequencies
    for (int c = 0; c < ncols; c++)
        for (int b = 0; b < MSA_ALPHA; b++)
            p->freq[c][b] /= (float)nseqs;
    return p;
}

// Score aligning two profile columns (vectors of frequencies).
// gap_penalty penalises gaps (index 5). The formula is:
//   score = sum_{b∈{A,C,G,T,N}} ca[b] * cb[b]  -  (ca[5] + cb[5]) * gap_penalty/2
static float profile_col_score(float *ca, float *cb, float gap_penalty) {
    float match = 0.0f;
    for (int b = 0; b < 5; b++) match += ca[b] * cb[b];
    return match - (ca[5] + cb[5]) * gap_penalty * 0.5f;
}

// ---------- Needleman-Wunsch profile-vs-profile ----------
#define NW_GAP (-1.5f)   // gap penalty for profile alignment

// A column pair in the alignment: ia = column index from first profile (or -1 for gap),
// ib = column index from second profile (or -1 for gap)
typedef struct { int ia; int ib; } ColPair;

// Perform Needleman‑Wunsch global alignment between two profiles.
// Returns an array of ColPair (length *out_len) describing the alignment in order.
static ColPair *nw_profile_align(Profile *pa, Profile *pb, int *out_len) {
    int la = pa->ncols, lb = pb->ncols;

    float *dp   = (float *)malloc((la+1)*(lb+1)*sizeof(float));
    int   *from = (int *)  malloc((la+1)*(lb+1)*sizeof(int));
    if (!dp || !from) { free(dp); free(from); return NULL; }

#define DP(i,j)   dp[(i)*(lb+1)+(j)]
#define FROM(i,j) from[(i)*(lb+1)+(j)]

    // Initialise DP matrix
    DP(0,0) = 0.0f;
    for (int i = 1; i <= la; i++) { DP(i,0) = NW_GAP * i; FROM(i,0) = 1; } // 1 = up (gap in pb)
    for (int j = 1; j <= lb; j++) { DP(0,j) = NW_GAP * j; FROM(0,j) = 2; } // 2 = left (gap in pa)

    // Fill DP
    for (int i = 1; i <= la; i++) {
        for (int j = 1; j <= lb; j++) {
            float s  = DP(i-1,j-1) + profile_col_score(pa->freq[i-1], pb->freq[j-1], -NW_GAP);
            float up = DP(i-1,j)   + NW_GAP;
            float lf = DP(i,  j-1) + NW_GAP;
            if      (s >= up && s >= lf) { DP(i,j) = s;  FROM(i,j) = 0; } // diagonal
            else if (up >= lf)           { DP(i,j) = up; FROM(i,j) = 1; } // up
            else                         { DP(i,j) = lf; FROM(i,j) = 2; } // left
        }
    }

    // Traceback to recover the alignment
    ColPair *pairs = (ColPair *)malloc((la + lb) * sizeof(ColPair));
    if (!pairs) { free(dp); free(from); return NULL; }

    int pi = 0, i = la, j = lb;
    while (i > 0 || j > 0) {
        int f = (i > 0 && j > 0) ? FROM(i,j) : (i > 0 ? 1 : 2);
        if      (f == 0) { pairs[pi++] = (ColPair){i-1, j-1}; i--; j--; }
        else if (f == 1) { pairs[pi++] = (ColPair){i-1,  -1}; i--;      } // gap in second profile
        else             { pairs[pi++] = (ColPair){ -1, j-1};      j--; } // gap in first profile
    }
    // Reverse to get correct order
    for (int a = 0, b = pi-1; a < b; a++, b--) {
        ColPair tmp = pairs[a]; pairs[a] = pairs[b]; pairs[b] = tmp;
    }

    free(dp); free(from);
    *out_len = pi;
    return pairs;

#undef DP
#undef FROM
}

// ---------- Aligned sequence set ----------
// Holds a collection of sequences that have been aligned to the same columns.
typedef struct {
    char **seqs;    // array of strings, each of length ncols (including gaps)
    int   nseqs;    // number of sequences
    int   ncols;    // aligned length (number of columns)
} AlignedSet;

static void aligned_set_free(AlignedSet *as) {
    if (!as) return;
    for (int i = 0; i < as->nseqs; i++) free(as->seqs[i]);
    free(as->seqs);
    free(as);
}

// Merge two aligned sets according to a column‑pair mapping (pairs).
// The resulting AlignedSet will have ncols = npairs, and each sequence is built
// by taking the appropriate column from the original set, inserting gaps when needed.
static AlignedSet *merge_aligned_sets(AlignedSet *a, AlignedSet *b,
                                      ColPair *pairs, int npairs) {
    AlignedSet *out = (AlignedSet *)malloc(sizeof(AlignedSet));
    if (!out) return NULL;
    out->nseqs = a->nseqs + b->nseqs;
    out->ncols = npairs;
    out->seqs  = (char **)malloc(out->nseqs * sizeof(char *));
    if (!out->seqs) { free(out); return NULL; }

    // Build each sequence column by column
    for (int s = 0; s < out->nseqs; s++) {
        out->seqs[s] = (char *)malloc(npairs + 1);
        if (!out->seqs[s]) {
            for (int k = 0; k < s; k++) free(out->seqs[k]);
            free(out->seqs); free(out); return NULL;
        }
        int from_a  = (s < a->nseqs);          // true for sequences originally from set a
        int si      = from_a ? s : (s - a->nseqs);
        char **src  = from_a ? a->seqs : b->seqs;
        for (int c = 0; c < npairs; c++) {
            int idx = from_a ? pairs[c].ia : pairs[c].ib;
            out->seqs[s][c] = (idx < 0) ? GAP_CHAR : src[si][idx];
        }
        out->seqs[s][npairs] = '\0';
    }
    return out;
}

// ---------- UPGMA guide tree ----------
// Each step records the indices of the two clusters being merged.
typedef struct { int ci; int cj; } GuideMerge;

// Build a UPGMA tree from a distance matrix (n x n, symmetric, diagonal 0).
// Returns an array of GuideMerge of length n-1 (or less on error).
static GuideMerge *upgma(double **dist, int n, int *out_nmerges) {
    GuideMerge *merges = (GuideMerge *)malloc((n-1) * sizeof(GuideMerge));
    double    **cd     = (double **)  malloc(n * sizeof(double *));
    int        *csz    = (int *)      calloc(n, sizeof(int));
    int        *alive  = (int *)      malloc(n * sizeof(int));
    if (!merges || !cd || !csz || !alive) {
        free(merges); free(cd); free(csz); free(alive); return NULL;
    }
    // Initialise: each cluster is one element, with its own copy of distances
    for (int i = 0; i < n; i++) {
        cd[i] = (double *)malloc(n * sizeof(double));
        if (!cd[i]) { free(merges); return NULL; }
        for (int j = 0; j < n; j++) cd[i][j] = dist[i][j];
        csz[i] = 1; alive[i] = 1;
    }

    int nm = 0;
    for (int step = 0; step < n-1; step++) {
        // Find the closest pair among alive clusters
        double best = DBL_MAX;
        int bi = -1, bj = -1;
        for (int i = 0; i < n; i++) {
            if (!alive[i]) continue;
            for (int j = i+1; j < n; j++) {
                if (!alive[j]) continue;
                if (cd[i][j] < best) { best = cd[i][j]; bi = i; bj = j; }
            }
        }
        if (bi < 0) break;
        merges[nm++] = (GuideMerge){bi, bj};

        // Update distances to the new merged cluster (stored in bi)
        int sa = csz[bi], sb = csz[bj];
        for (int k = 0; k < n; k++) {
            if (!alive[k] || k == bi || k == bj) continue;
            double d = (cd[bi][k]*sa + cd[bj][k]*sb) / (double)(sa+sb);
            cd[bi][k] = cd[k][bi] = d;
        }
        csz[bi] = sa + sb;
        alive[bj] = 0;  // mark bj as merged (dead)
    }

    for (int i = 0; i < n; i++) free(cd[i]);
    free(cd); free(csz); free(alive);
    *out_nmerges = nm;
    return merges;
}

// ---------- Consensus from a finished alignment ----------
// Build a consensus string from an aligned set:
//   - Columns where gaps dominate (>50%) are removed entirely.
//   - For remaining columns, the most frequent base (A/C/G/T/N) is chosen.
static char *consensus_from_alignment(AlignedSet *as) {
    if (!as || as->ncols == 0 || as->nseqs == 0) return NULL;
    char *cons = (char *)malloc(as->ncols + 1);
    if (!cons) return NULL;

    int ci = 0;
    for (int c = 0; c < as->ncols; c++) {
        int counts[MSA_ALPHA] = {0};
        for (int s = 0; s < as->nseqs; s++)
            counts[base_index(as->seqs[s][c])]++;

        // Skip columns where more than half of the sequences have a gap
        if (counts[5] * 2 > as->nseqs) continue;

        int best = 0, bi = 4;  // bi = 4 (N) is fallback
        for (int b = 0; b < 5; b++)   // only consider A,C,G,T,N (ignore gap)
            if (counts[b] > best) { best = counts[b]; bi = b; }
        cons[ci++] = index_base(bi);
    }
    cons[ci] = '\0';
    return cons;
}

// ---------- MSA for one connected component ----------
// seqs[] : array of raw (unaligned) sequences in the component
// n      : number of sequences
// merges_arr[] : all merge possibilities discovered earlier (global list)
// nmerges_arr   : number of entries in merges_arr
// seq_indices[] : mapping from local index (0..n-1) to global contig index
// Returns a consensus string (heap-allocated) for this component.
static char *msa_consensus(char **seqs, int n,
                            merge_info *merges_arr, int nmerges_arr,
                            int *seq_indices) {
    if (n == 0) return NULL;
    if (n == 1) return strdup(seqs[0]);

    // 1. Build pairwise distance matrix (n x n) from the substitution rates observed in overlaps.
    double **dist = (double **)malloc(n * sizeof(double *));
    if (!dist) return NULL;
    for (int i = 0; i < n; i++) {
        dist[i] = (double *)malloc(n * sizeof(double));
        if (!dist[i]) { for (int k=0;k<i;k++) free(dist[k]); free(dist); return NULL; }
        for (int j = 0; j < n; j++) dist[i][j] = (i == j) ? 0.0 : 1.0; // default = 1 (max distance)
    }
    // Fill distances using the merge_info records that involve these sequences
    for (int mi = 0; mi < nmerges_arr; mi++) {
        merge_info info = merges_arr[mi];
        if (info.merged_length == 0) continue;
        int pi = -1, si = -1;
        for (int k = 0; k < n; k++) {
            if (seq_indices[k] == info.index_p) pi = k;
            if (seq_indices[k] == info.index_s) si = k;
        }
        if (pi < 0 || si < 0) continue;
        // distance = mismatches / overlap length
        double d = (double)info.substitutions /
                   (double)(info.overlap_length > 0 ? info.overlap_length : 1);
        if (d < dist[pi][si]) dist[pi][si] = dist[si][pi] = d;
    }

    // 2. UPGMA guide tree from the distance matrix
    int ngt = 0;
    GuideMerge *guide = upgma(dist, n, &ngt);
    for (int i = 0; i < n; i++) free(dist[i]);
    free(dist);
    if (!guide) return NULL;

    // 3. Initialise each sequence as its own aligned set (trivial alignment)
    AlignedSet **clusters = (AlignedSet **)malloc(n * sizeof(AlignedSet *));
    if (!clusters) { free(guide); return NULL; }
    for (int i = 0; i < n; i++) {
        clusters[i] = (AlignedSet *)malloc(sizeof(AlignedSet));
        if (!clusters[i]) { free(guide); free(clusters); return NULL; }
        clusters[i]->nseqs   = 1;
        clusters[i]->ncols   = (int)strlen(seqs[i]);
        clusters[i]->seqs    = (char **)malloc(sizeof(char *));
        clusters[i]->seqs[0] = strdup(seqs[i]);
    }

    // 4. Progressive alignment following the guide tree (bottom-up)
    for (int step = 0; step < ngt; step++) {
        int ci = guide[step].ci, cj = guide[step].cj;
        if (!clusters[ci] || !clusters[cj]) continue;

        // Build profiles for the two clusters
        Profile *pa = profile_from_seqs(clusters[ci]->seqs, clusters[ci]->nseqs, clusters[ci]->ncols);
        Profile *pb = profile_from_seqs(clusters[cj]->seqs, clusters[cj]->nseqs, clusters[cj]->ncols);
        if (!pa || !pb) { profile_free(pa); profile_free(pb); continue; }

        // Align the two profiles
        int npairs = 0;
        ColPair *pairs = nw_profile_align(pa, pb, &npairs);
        profile_free(pa); profile_free(pb);
        if (!pairs) continue;

        // Merge the two aligned sets using the obtained column mapping
        AlignedSet *merged = merge_aligned_sets(clusters[ci], clusters[cj], pairs, npairs);
        free(pairs);
        if (!merged) continue;

        // Replace the two old clusters with the new merged one
        aligned_set_free(clusters[ci]);
        aligned_set_free(clusters[cj]);
        clusters[ci] = merged;
        clusters[cj] = NULL;
    }

    // 5. Extract consensus from the root cluster (the last non‑NULL one)
    AlignedSet *root = NULL;
    for (int i = 0; i < n; i++) if (clusters[i]) { root = clusters[i]; break; }
    char *cons = consensus_from_alignment(root);

    // Cleanup
    for (int i = 0; i < n; i++) if (clusters[i]) aligned_set_free(clusters[i]);
    free(clusters);
    free(guide);
    return cons;
}

// =============================================================================
// filter_contained_contigs (parallelized)
// =============================================================================
// Remove sequences that are near-duplicates / substrings of another (longer) sequence
// in the list, allowing up to f_max_substitutions_dedup mismatches. Contigs shorter
// than min_length are skipped entirely (neither compared as a candidate nor removed).
// Returns a new array of unique, non‑contained sequences.

// ---------- Thread data for parallel containment filtering ----------
typedef struct {
    int thread_id;
    int number_threads;
    int count;
    char **contigs;
    int *lengths;
    int *is_contained;
} contain_thread_data_t;

static void *filter_contained_worker(void *arg) {
    contain_thread_data_t *data = (contain_thread_data_t *)arg;
    int thread_id     = data->thread_id;
    int nt            = data->number_threads;
    int count         = data->count;
    char **contigs    = data->contigs;
    int *lengths      = data->lengths;
    int *is_contained = data->is_contained;

    // Partition unordered pairs (i < j) across threads by i, so each pair is
    // examined exactly once instead of twice (i,j) and (j,i).
    for (int i = thread_id; i < count; i += nt) {
        int li = lengths[i];
        if (li < min_length) continue;

        for (int j = i + 1; j < count; j++) {
            int lj = lengths[j];
            if (lj < min_length) continue;

            if (li == lj) {
                int max_subs = (int)(f_max_substitutions_dedup * li);
                if (hamming_distance_limit(contigs[i], contigs[j], li, max_subs) <= max_subs) {
                    // Same length + near-identical: mark the later index as the duplicate.
                    __atomic_store_n(&is_contained[j], 1, __ATOMIC_RELAXED);
                }
                continue;
            }

            const char *longer, *shorter;
            int shorter_len, *mark_idx;
            if (li > lj) { longer = contigs[i]; shorter = contigs[j]; shorter_len = lj; mark_idx = &is_contained[j]; }
            else         { longer = contigs[j]; shorter = contigs[i]; shorter_len = li; mark_idx = &is_contained[i]; }

            int max_subs = (int)(f_max_substitutions_dedup * shorter_len);
            if (contains_with_substitutions(longer, shorter, max_subs)) {
                __atomic_store_n(mark_idx, 1, __ATOMIC_RELAXED);
            }
        }
    }
    return NULL;
}

char **filter_contained_contigs(char **contigs, int count, int *new_count) {
    if (count == 0) { *new_count = 0; return NULL; }

    int *is_contained = (int *)calloc(count, sizeof(int));
    int *lengths       = (int *)malloc(count * sizeof(int));
    if (!is_contained || !lengths) { free(is_contained); free(lengths); return NULL; }
    for (int i = 0; i < count; i++) lengths[i] = (int)strlen(contigs[i]); // cache once

    int nt = number_threads > 0 ? number_threads : 1;
    if (nt > count) nt = count;
    if (nt < 1) nt = 1;

    pthread_t threads[nt];
    contain_thread_data_t tdata[nt];
    for (int t = 0; t < nt; t++) {
        tdata[t] = (contain_thread_data_t){
            .thread_id = t, .number_threads = nt, .count = count,
            .contigs = contigs, .lengths = lengths, .is_contained = is_contained
        };
        pthread_create(&threads[t], NULL, filter_contained_worker, &tdata[t]);
    }
    for (int t = 0; t < nt; t++) pthread_join(threads[t], NULL);

    free(lengths);

    *new_count = 0;
    for (int i = 0; i < count; i++) if (!is_contained[i]) (*new_count)++;

    char **filtered = (char **)malloc((*new_count) * sizeof(char *));
    if (!filtered) { free(is_contained); return NULL; }

    int idx = 0;
    for (int i = 0; i < count; i++) {
        if (!is_contained[i]) {
            filtered[idx] = strdup(contigs[i]);
            if (!filtered[idx]) {
                for (int k = 0; k < idx; k++) free(filtered[k]);
                free(filtered); free(is_contained); return NULL;
            }
            idx++;
        }
    }
    free(is_contained);
    return filtered;
}

// =============================================================================
// Parallel per-component MSA worker (used in merge_contigs Phase 3)
// =============================================================================
typedef struct {
    int thread_id;
    int number_threads;
    int nroots;
    int n;
    int *roots;
    int *comp_of;
    merge_info *merges_arr;
    int nmerges_arr;
    char **consensi;   // output slots, one per root index r (NULL if none produced)
} msa_thread_data_t;

static void *msa_worker(void *arg) {
    msa_thread_data_t *data = (msa_thread_data_t *)arg;
    int thread_id = data->thread_id;
    int nt        = data->number_threads;

    for (int r = thread_id; r < data->nroots; r += nt) {
        int root = data->roots[r];

        // Collect indices of all contigs belonging to this component
        int mem_cap = 8, mem_cnt = 0;
        int *members = (int *)malloc(mem_cap * sizeof(int));
        if (!members) continue;
        for (int i = 0; i < data->n; i++) {
            if (data->comp_of[i] != root) continue;
            if (mem_cnt == mem_cap) {
                mem_cap *= 2;
                int *tmp = (int *)realloc(members, mem_cap * sizeof(int));
                if (!tmp) break;
                members = tmp;
            }
            members[mem_cnt++] = i;
        }
        if (mem_cnt == 0) { free(members); continue; }

        // Build an array of sequence strings for this component
        char **comp_seqs = (char **)malloc(mem_cnt * sizeof(char *));
        if (!comp_seqs) { free(members); continue; }
        for (int k = 0; k < mem_cnt; k++) comp_seqs[k] = list_contigs[members[k]].contig;

        if (verbose)
            printf("Thread %d: component root=%d: %d sequences -> running MSA\n",
                   thread_id, root, mem_cnt);

        // Compute consensus for this component using MSA (each thread only ever
        // touches its own comp_seqs/members and writes to its own consensi[r] slot,
        // so no locking is required here)
        char *cons = msa_consensus(comp_seqs, mem_cnt, data->merges_arr, data->nmerges_arr, members);
        free(comp_seqs);
        free(members);

        if (cons && strlen(cons) > 0) {
            data->consensi[r] = cons;
        } else {
            data->consensi[r] = NULL;
            free(cons);
        }
    }
    return NULL;
}

// =============================================================================
// merge_contigs
// =============================================================================
// Main entry point: given the global list_contigs (size number_sequences),
// finds all possible merges (allowing mismatches), groups contigs into connected components,
// builds a consensus for each component via progressive MSA, and finally filters out
// contained sequences. Returns a contig_array with the merged contigs.
contig_array merge_contigs(int number_threads, int number_sequences) {

    printf("Starting merge of contigs with %d threads...\n", number_threads);

    if (number_sequences <= 0) {
        printf("Invalid number of sequences.\n");
        return (contig_array){NULL, 0};
    }

    long max_pairs = ((long)number_sequences * (number_sequences - 1)) / 2;
    if (max_pairs < 1) max_pairs = 1;

    merge_possibilities = (merge_info *)malloc(max_pairs * sizeof(merge_info));
    if (!merge_possibilities) {
        printf("Memory allocation failed for merge possibilities.\n");
        return (contig_array){NULL, 0};
    }
    merge_capacity              = (int)max_pairs;
    index_merge_possibilities   = 0;

    // ---- Phase 1: gather pairwise overlaps in parallel ----
    pthread_t     threads[number_threads];
    thread_data_t thread_data[number_threads];
    for (int t = 0; t < number_threads; t++) {
        thread_data[t].thread_id       = t;
        thread_data[t].number_threads  = number_threads;
        thread_data[t].total_sequences = number_sequences;
        pthread_create(&threads[t], NULL, get_merge_possibilities, &thread_data[t]);
    }
    for (int t = 0; t < number_threads; t++) pthread_join(threads[t], NULL);

    int nmp = index_merge_possibilities;
    if (nmp == 0) {
        printf("No merge possibilities found.\n");

        // Allocate an array of pointers (copy of list_contigs)
        char **contigs_out = malloc(number_sequences * sizeof(char *));
        if (contigs_out == NULL) {
            // Allocation failed; return an empty result (or handle as needed)
            return (contig_array){NULL, 0};
        }

        for (int i = 0; i < number_sequences; i++) {
            contigs_out[i] = list_contigs[i].contig;
        }

        int filtered_count = 0;
        char **filtered = filter_contained_contigs(contigs_out, number_sequences, &filtered_count);

        // Free the temporary copy; filtered contains the result (assumed newly allocated)
        free(contigs_out);

        return (contig_array){filtered, filtered_count};
    }

    // ---- Phase 2: Union-Find to identify connected components ----
    int  n       = number_sequences;
    int *parent  = (int *)malloc(n * sizeof(int));
    int *rank_uf = (int *)calloc(n, sizeof(int));
    if (!parent || !rank_uf) {
        free(parent); free(rank_uf);
        free(merge_possibilities); merge_possibilities = NULL;
        return (contig_array){NULL, 0};
    }
    for (int i = 0; i < n; i++) parent[i] = i;

#define UF_FIND(x) ({ int _x=(x); while(parent[_x]!=_x){parent[_x]=parent[parent[_x]];_x=parent[_x];} _x; })

    // Union the two contigs for each valid merge possibility
    for (int mi = 0; mi < nmp; mi++) {
        merge_info info = merge_possibilities[mi];
        if (info.merged_length == 0) continue;
        int rp = UF_FIND(info.index_p), rs = UF_FIND(info.index_s);
        if (rp != rs) {
            if (rank_uf[rp] < rank_uf[rs]) { int tmp=rp; rp=rs; rs=tmp; }
            parent[rs] = rp;
            if (rank_uf[rp] == rank_uf[rs]) rank_uf[rp]++;
        }
    }

    // Find the root for each contig
    int *comp_of = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) comp_of[i] = UF_FIND(i);

    // Collect unique roots
    int *roots  = (int *)malloc(n * sizeof(int));
    int  nroots = 0;
    for (int i = 0; i < n; i++) {
        int found = 0;
        for (int k = 0; k < nroots; k++) if (roots[k] == comp_of[i]) { found=1; break; }
        if (!found) roots[nroots++] = comp_of[i];
    }
    if (verbose) printf("Connected components: %d\n", nroots);

    // ---- Phase 3: MSA per component -> consensus (parallelized across components) ----
    char **consensi_slots = (char **)calloc(nroots, sizeof(char *)); // one slot per root, NULL if empty
    if (!consensi_slots) {
        free(parent); free(rank_uf); free(comp_of); free(roots);
        free(merge_possibilities); merge_possibilities = NULL;
        return (contig_array){NULL, 0};
    }

    int nt_msa = number_threads > 0 ? number_threads : 1;
    if (nt_msa > nroots) nt_msa = nroots > 0 ? nroots : 1;
    if (nt_msa < 1) nt_msa = 1;

    pthread_t msa_threads[nt_msa];
    msa_thread_data_t msa_tdata[nt_msa];
    for (int t = 0; t < nt_msa; t++) {
        msa_tdata[t] = (msa_thread_data_t){
            .thread_id = t, .number_threads = nt_msa, .nroots = nroots, .n = n,
            .roots = roots, .comp_of = comp_of,
            .merges_arr = merge_possibilities, .nmerges_arr = nmp,
            .consensi = consensi_slots
        };
        pthread_create(&msa_threads[t], NULL, msa_worker, &msa_tdata[t]);
    }
    for (int t = 0; t < nt_msa; t++) pthread_join(msa_threads[t], NULL);

    char **consensi  = (char **)malloc(nroots * sizeof(char *));
    int   cons_count = 0;
    if (!consensi) {
        for (int r = 0; r < nroots; r++) free(consensi_slots[r]);
        free(consensi_slots);
        free(parent); free(rank_uf); free(comp_of); free(roots);
        free(merge_possibilities); merge_possibilities = NULL;
        return (contig_array){NULL, 0};
    }
    for (int r = 0; r < nroots; r++) {
        if (consensi_slots[r]) consensi[cons_count++] = consensi_slots[r];
    }
    free(consensi_slots);

    free(parent); free(rank_uf); free(comp_of); free(roots);
    free(merge_possibilities); merge_possibilities = NULL;

    // ---- Phase 4: filter contained contigs (remove substrings, parallelized) ----
    int   filtered_count = 0;
    char **filtered = filter_contained_contigs(consensi, cons_count, &filtered_count);
    for (int i = 0; i < cons_count; i++) free(consensi[i]);
    free(consensi);

    printf("-----Total contigs after merging: %d\n", filtered_count);
    return (contig_array){filtered, filtered_count};
}

// ---------- Try to merge two contigs (optimised) ----------
// Same as try_merge_two_contigs, but actually builds and returns the merged string.
// Returns a new allocated string on success, NULL otherwise.
char *try_merge_two_reads(char *a, int len_a, char *b, int len_b, float perc_max_subs, int min_overlap) {
    int max_possible = len_a < len_b ? len_a : len_b;
    if (min_overlap > max_possible) return NULL;

    int overlap_start = max_possible;
    if (max_overlap > 0 && overlap_start > max_overlap)
        overlap_start = max_overlap;

    for (int overlap = overlap_start; overlap >= min_overlap; --overlap) {

        int max_subs = perc_max_subs * avg_read_length;
        
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