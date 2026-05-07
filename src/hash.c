#include <stdio.h>
#include <stdlib.h>
#include <string.h>      // for strdup, strcmp
#include "mem.h"         // provides Calloc, Realloc, Free
#include "hash.h"

// ---------- Hash function for strings (djb2 + ZHASH mixing) ----------
static uint64_t ZHASH(uint64_t z) {
    z = (~z) + (z << 21);
    z = z    ^ (z >> 24);
    z = (z   + (z << 3)) + (z << 8);
    z = z    ^ (z >> 14);
    z = (z   + (z << 2)) + (z << 4);
    z = z    ^ (z >> 28);
    z = z    + (z << 31);
    return z;
}

static uint64_t HashString(const char *s) {
    uint64_t h = 5381;
    int c;
    while ((c = *s++))
        h = ((h << 5) + h) + c;
    return ZHASH(h);
}

// ---------- Create an empty hash table ----------
HASH *CreateHashTable(void) {
    HASH *H = (HASH *) Calloc(1, sizeof(HASH));
    H->size      = HASH_SIZE;
    H->entries   = (ENTRY  **) Calloc(H->size, sizeof(ENTRY *));
    H->entrySize = (uint64_t *) Calloc(H->size, sizeof(uint64_t));
    return H;
}

// ---------- Insert a new entry into a specific bucket ----------
void InsertEntry(HASH *H, uint64_t hi, const char *key, int value) {
    // Expand the bucket array
    H->entries[hi] = (ENTRY *) Realloc(H->entries[hi],
                                       (H->entrySize[hi] + 1) * sizeof(ENTRY));
    // Create the new entry
    ENTRY *newEntry = &H->entries[hi][H->entrySize[hi]];
    newEntry->key   = strdup(key);      // duplicate the string
    newEntry->value = value;
    H->entrySize[hi]++;
}

// ---------- Retrieve the value for a given key ----------
// Returns a pointer to the value if found, otherwise NULL.
int *GetValue(HASH *H, const char *key) {
    uint64_t hash = HashString(key);
    uint64_t hi   = hash % H->size;

    for (uint64_t x = 0; x < H->entrySize[hi]; ++x) {
        if (strcmp(H->entries[hi][x].key, key) == 0)
            return &H->entries[hi][x].value;
    }
    return NULL;
}

// ---------- Insert or update a key-value pair ----------
void SetValue(HASH *H, const char *key, int value) {
    uint64_t hash = HashString(key);
    uint64_t hi   = hash % H->size;

    // Search for existing key
    for (uint64_t x = 0; x < H->entrySize[hi]; ++x) {
        if (strcmp(H->entries[hi][x].key, key) == 0) {
            H->entries[hi][x].value = value;   // update
            return;
        }
    }
    // Key not found → insert new entry
    InsertEntry(H, hi, key, value);
}

// ---------- Free the entire hash table ----------
void RemoveHashTable(HASH *H) {
    for (uint64_t i = 0; i < H->size; ++i) {
        if (H->entrySize[i] != 0) {
            for (uint64_t j = 0; j < H->entrySize[i]; ++j)
                Free(H->entries[i][j].key);    // free the duplicated string
            Free(H->entries[i]);               // free the bucket array
        }
    }
    Free(H->entries);
    Free(H->entrySize);
    Free(H);
}