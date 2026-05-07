#ifndef HASH_H
#define HASH_H

#include <stdint.h>

#define HASH_SIZE 1024

typedef struct entry {
    char *key;
    int  value;
} ENTRY;

typedef struct hash {
    uint64_t size;
    ENTRY   **entries;
    uint64_t *entrySize;
} HASH;

HASH  *CreateHashTable(void);
void   SetValue(HASH *H, const char *key, int value);
int   *GetValue(HASH *H, const char *key);
void   RemoveHashTable(HASH *H);

#endif