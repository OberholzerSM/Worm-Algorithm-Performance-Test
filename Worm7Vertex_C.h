#pragma once

//If including the header from the main C++ file, the functions need to be marked as having been compiled in C.
#ifdef __cplusplus
extern "C" {
#else
//Have only the C file include the libraries and define the BondType enum.
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
//BondTypes
enum BondType
{
	NO_BOND,
	BOND,
};

//PCG Struct
struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
};
typedef struct pcg_state_setseq_64 pcg32_random_t;

//Inital Seed
#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

void pcg32_srandom(uint64_t initstate, uint64_t initseq);
void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate,
                     uint64_t initseq);

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

uint32_t pcg32_random(void);
uint32_t pcg32_random_r(pcg32_random_t* rng);

// pcg32_boundedrand(bound):
// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound

uint32_t pcg32_boundedrand(uint32_t bound);
uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound);

#endif

//Generates n_data many closed loop configurations.
void generateWormData_7Vertex_C(const int N0, const int N1, const double M, const int n_data);

//Generates a single closed loop configuration.
void Worm_7Vertex_C(const int N0, const int N1, const double M, enum BondType* const lattice);

//Determines the weight at a vertex. bool sink: If a sink/source is present.
double weight_7Vertex_C(const double M, const enum BondType vertex[4], const bool sink);

//Determines the proposal probability of moving the worm from x -> y.
double proposal_7Vertex_C(const int x[2], const int y[2], const enum BondType oldvertex_x[4], const double p0, const int N0, const int N1);

//Return a random int between and including a and b.
inline int randomint_C(const int a, const int b);

//Return a random double between and including 0 and 1.
inline double randomdouble_C();

#ifdef __cplusplus
}
#endif