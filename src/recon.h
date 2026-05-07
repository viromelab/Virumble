#ifndef CM_HASH_H_INCLUDED
#define CM_HASH_H_INCLUDED

#include "defs.h"

typedef uint8_t  HCC;        // Size of context counters for hash tables
typedef uint16_t ENTMAX;     // Entry size / capacity per hash index

#define MAX_HASH_COUNTER ((1<<(sizeof(HCC)*8))-1)
#define HASH_SIZE 16777259

typedef struct{
  int used;                  // Whether this entry is used
  int init_pos;              // Initial position of the sequence
  int end_pos;               // End position of the sequence
  }
SEQ_INFO;

typedef struct{
  
  }
HASH;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static inline uint64_t ZHASH(uint64_t z)
  {
  z = (~z) + (z << 21);
  z = z ^ (z >> 24);
  z = (z + (z << 3)) + (z << 8);
  z = z ^ (z >> 14);
  z = (z + (z << 2)) + (z << 4);
  z = z ^ (z >> 28);
  z = z + (z << 31);
  return z;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HASH    *CreateHashTable   (uint32_t);
void    InsertKey          (HASH *, uint32_t, uint64_t, uint8_t);
HCC     *GetHCCounters     (HASH *, uint64_t);
void    UpdateHashCounter  (HASH *, uint32_t, uint64_t);
void    RemoveHashTable    (HASH *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif