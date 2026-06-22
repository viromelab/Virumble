#ifndef AUX_H
#define AUX_H

#include <stdio.h>

/* Enum for number type */
typedef enum {
    TYPE_INT,
    TYPE_FLOAT,
    TYPE_INVALID
} NumType;

/**
 * Checks whether two paired‑end FASTQ files are synchronized.
 *
 * @param file1       Path to the first FASTQ file (read 1)
 * @param file2       Path to the second FASTQ file (read 2)
 * @param error_msg   Buffer to store an error message (can be NULL)
 * @param msg_size    Size of error_msg buffer (ignored if error_msg is NULL)
 *
 * @return 0 if files are synchronized,
 *         1 if file count mismatch or identifier mismatch,
 *         2 if a file cannot be opened,
 *         3 if a record is incomplete (truncated file).
 */
int check_paired_end_files(const char *file1, const char *file2,
                           char *error_msg, size_t msg_size);
NumType parse_number(const char *str, double *value);

#endif /* AUX_H */