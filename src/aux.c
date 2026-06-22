
#include "aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#define MAX_LINE 4096

/* Normalize a FASTQ identifier:
   - remove leading '@'
   - truncate at first space/tab
   - remove trailing "/1" or "/2"
   Result is stored in 'out' (must be at least MAX_LINE bytes). */
static void normalize_id(const char *raw, char *out) {
    const char *p = raw;
    char *dst = out;

    if (*p == '@')
        p++;

    while (*p && *p != ' ' && *p != '\t' && *p != '\n')
        *dst++ = *p++;
    *dst = '\0';

    size_t len = strlen(out);
    if (len >= 2 && out[len-2] == '/' && (out[len-1] == '1' || out[len-1] == '2'))
        out[len-2] = '\0';
}

/* Read one FASTQ record from fp.
   Returns: 1 = success, 0 = EOF, -1 = incomplete record. */
static int read_fastq_record(FILE *fp, char *id_buf) {
    char line[MAX_LINE];

    if (fgets(line, sizeof(line), fp) == NULL)
        return 0;
    line[strcspn(line, "\n")] = '\0';
    normalize_id(line, id_buf);

    if (fgets(line, sizeof(line), fp) == NULL) return -1;
    if (fgets(line, sizeof(line), fp) == NULL) return -1;
    if (fgets(line, sizeof(line), fp) == NULL) return -1;

    return 1;
}

int check_paired_end_files(const char *file1, const char *file2,
                           char *error_msg, size_t msg_size) {
    FILE *fp1 = fopen(file1, "r");
    if (!fp1) {
        if (error_msg)
            snprintf(error_msg, msg_size, "Cannot open file: %s", file1);
        return 2;
    }
    FILE *fp2 = fopen(file2, "r");
    if (!fp2) {
        fclose(fp1);
        if (error_msg)
            snprintf(error_msg, msg_size, "Cannot open file: %s", file2);
        return 2;
    }

    char id1[MAX_LINE], id2[MAX_LINE];
    long rec_num = 1;
    int ret1, ret2;

    while (1) {
        ret1 = read_fastq_record(fp1, id1);
        ret2 = read_fastq_record(fp2, id2);

        if (ret1 == 0 && ret2 == 0)   /* both EOF -> OK */
            break;

        if (ret1 < 0 || ret2 < 0) {
            fclose(fp1); fclose(fp2);
            if (error_msg)
                snprintf(error_msg, msg_size,
                         "Incomplete record at record %ld", rec_num);
            return 3;
        }

        if (ret1 == 0) {
            fclose(fp1); fclose(fp2);
            if (error_msg)
                snprintf(error_msg, msg_size,
                         "File2 has more reads than file1 (mismatch at record %ld)", rec_num);
            return 1;
        }
        if (ret2 == 0) {
            fclose(fp1); fclose(fp2);
            if (error_msg)
                snprintf(error_msg, msg_size,
                         "File1 has more reads than file2 (mismatch at record %ld)", rec_num);
            return 1;
        }

        if (strcmp(id1, id2) != 0) {
            fclose(fp1); fclose(fp2);
            if (error_msg)
                snprintf(error_msg, msg_size,
                         "Identifier mismatch at record %ld: '%s' vs '%s'",
                         rec_num, id1, id2);
            return 1;
        }

        rec_num++;
    }

    fclose(fp1);
    fclose(fp2);
    if (error_msg)
        snprintf(error_msg, msg_size, "Synchronized (checked %ld records)", rec_num - 1);
    return 0;
}


/* Trim whitespace in-place */
static void trim(char *str) {
    char *end;
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return;
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    *(end + 1) = '\0';
}

/*
 * Parse a string and determine if it's an integer or a float.
 * Stores the numeric value in *value (if valid).
 * Returns TYPE_INT, TYPE_FLOAT, or TYPE_INVALID.
 */
NumType parse_number(const char *str, double *value) {
    char *endptr;
    double val;
    char *tmp = NULL;
    int has_decimal = 0, has_exponent = 0;
    const char *p;

    /* Work on a copy so we can trim without modifying the original */
    tmp = strdup(str);
    if (!tmp) return TYPE_INVALID;
    trim(tmp);

    if (*tmp == '\0') {
        free(tmp);
        return TYPE_INVALID;
    }

    /* Detect presence of '.' or 'e'/'E' */
    p = tmp;
    if (*p == '+' || *p == '-') p++;
    while (*p) {
        if (*p == '.') has_decimal = 1;
        else if (*p == 'e' || *p == 'E') has_exponent = 1;
        p++;
    }

    errno = 0;
    val = strtod(tmp, &endptr);
    if (errno == ERANGE || *endptr != '\0') {
        free(tmp);
        return TYPE_INVALID;
    }

    free(tmp);

    if (has_decimal || has_exponent) {
        *value = val;
        return TYPE_FLOAT;
    } else {
        /* It's a plain integer (possibly with sign) */
        *value = val;
        return TYPE_INT;
    }
}