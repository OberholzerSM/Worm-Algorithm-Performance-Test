#include "Worm7Vertex_CUDA.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//PCG Struct
struct pcg_state_setseq_64 
{
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
};
typedef struct pcg_state_setseq_64 pcg32_random_t;

//Inital Seed
#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }

//PCG, code taken from https://github.com/imneme/pcg-c-basic (23.02.2026)

// state for global RNGs
static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

__device__ uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = (uint32_t)( ((oldstate >> 18u) ^ oldstate) >> 27u );
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((0-rot) & 31));
}

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

__device__ void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

// pcg32_boundedrand(bound):
// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound

__device__ uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)
{
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_t threshold = (0-bound) % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In 
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;) {
        uint32_t r = pcg32_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}

__device__ int randomint_Cuda(const int a, const int b, pcg32_random_t* rng)
{
	const uint32_t d = (uint32_t)(b-a+1);
	const uint32_t r = pcg32_boundedrand_r(rng,d);
	return (int)r + a;
}

__device__ double randomdouble_Cuda(pcg32_random_t* rng)
{
	return (double)pcg32_random_r(rng) / (double)UINT32_MAX;
}

//Worm Functions

enum BondType
{
	NO_BOND,
	BOND,
};

//Macro to more easily access array elements.
#define LATTICE(i,j,k) lattice[k + 4*j + 4*N1*i]

//Constant
#define ROOT3_INV4 0.6299605249474366

__device__ double weight_7Vertex_Cuda(double M, const enum BondType vertex[4], bool sink)
{
	//Determine number of bonds at vertex.
	int n_bonds = 0;
	for(int k=0;k<4;k++)
	{
		if(vertex[k] != NO_BOND)
			++n_bonds;
	}

	switch(n_bonds)
	{
	case 0:
		return ( sink ? 4.0 : M*M );

	case 1:
		return ( sink ? 1.0 : 0.0 );

	case 2:
	{
		const bool straightLine = ( vertex[0] == BOND && vertex[2] == BOND ) || ( vertex[1] == BOND && vertex[3] == BOND );
		return ( straightLine ? 1.0 : 0.5 );
	}

	case 3:
		return ( sink ? ROOT3_INV4 : 0.0 );

	default:
		return 0.0;
	}
}

__device__ double proposal_7Vertex_Cuda(const int x[2], const int y[2], const enum BondType oldvertex_x[4], const double p0, const int N0, const int N1)
{
	if(x[0]==N0)
		return 1.0 / (double)(N0*N1);
	else
	{
		//Determine number of vonds at oldvertex_x.
		int n_bonds = 0;
		for(int k=0;k<4;k++)
		{
			if(oldvertex_x[k] != NO_BOND)
				++n_bonds;
		}

		switch(n_bonds)
		{
		//Wormhead can be removed.
		case 0:
		case 2:
			return ( (y[0]==N0) ? p0 : (1.0-p0) / 4.0 );

		//Wormhead cannot be removed.
		case 1:
		case 3:
			return 0.25;

		default:
			return 0.0;
		}
	}

	return 0.0;
}

__device__ void Worm_7Vertex_Cuda(const int N0, const int N1, const double M, enum BondType* const lattice, pcg32_random_t* rng)
{
	//Suggest a new starting location for the worm.
	const int start[2] = { randomint_Cuda(0,N0-1,rng), randomint_Cuda(0,N1-1,rng) };
	int x[2];

	//Probablity to propose removing the worm.
	//const double p0 = 1.0 / (double)(N0*N1);
	const double p0 = 0.5;

	//See if the start will be accepted.
	bool check_end = true;
	{
		//Determine the acceptance probability.
		enum BondType vertex_start[4];
		for(int k=0; k<4; k++)
			vertex_start[k] = LATTICE(start[0],start[1],k);
		const double weightx = weight_7Vertex_Cuda(M,vertex_start,false);
		const double weightx_new = weight_7Vertex_Cuda(M,vertex_start,true);
		const double q1 = 1.0 / (double)(N0*N1);
		const double q2 = p0;
		const double pA = (weightx_new / weightx) * (q2/q1);

		//Test if the start will be accepted.
		const double u = randomdouble_Cuda(rng);
		if(u < pA)
		{
			check_end = false;
			x[0] = start[0];
			x[1] = start[1];
		}
	}

	//If the start was accepted, begin the worm loop.
	while(!check_end)
	{
		//Set the current vertex at location x.
		enum BondType oldvertex_x[4];
		for(int k=0; k<4; k++)
			oldvertex_x[k] = LATTICE(x[0],x[1],k);
		enum BondType newvertex_x[4];

		//New location y.
		int y[2];

		//Vertices at y.
		enum BondType oldvertex_y[4];
		enum BondType newvertex_y[4];

		//If a loop has been completed, propose to remove the worm.
		bool check_end_prop = false;
		if( x[0] == start[0] && x[1] == start[1] )
		{
			const double u = randomdouble_Cuda(rng);
			if(u < p0)
			{
				check_end_prop = true;
				y[0] = N0;
				y[1] = N1;
				for(int k=0; k<4; k++)
				{
					newvertex_x[k] = oldvertex_x[k];
					oldvertex_y[k] = oldvertex_x[k];
					newvertex_y[k] = oldvertex_x[k];
				}
			}
		}

		//Otherwise, propose a new location for the worm to move to.
		if(!check_end_prop)
		{
			//Propose a new direction to move towards.
			enum Direction
			{
				RIGHT,
				UP,
				LEFT,
				DOWN,
			};
			const enum Direction direction = (enum Direction)randomint_Cuda(0,3,rng);

			switch(direction)
			{
			case(RIGHT):

				//Define the new location y.
				if( x[0] < N0-1 )
					y[0] = x[0] + 1;
				else
					y[0] = 0;
				y[1] = x[1];

				//Define vertices.
				for(int k=0; k<4; k++)
				{
					oldvertex_y[k] = LATTICE(y[0],y[1],k);
					newvertex_x[k] = oldvertex_x[k];
					newvertex_y[k] = oldvertex_y[k];
				}

				//If no bond exists, create a new one.
				if( oldvertex_x[RIGHT] == NO_BOND )
				{
					newvertex_x[RIGHT] = BOND;
					newvertex_y[LEFT] = BOND;
				}
				//If a bond already exists, delete it.
				else
				{
					newvertex_x[RIGHT] = NO_BOND;
					newvertex_y[LEFT] = NO_BOND;
				}

				break;

			case(UP):

				//Define the new location y.
				if( x[1] < N1-1 )
					y[1] = x[1] + 1;
				else
					y[1] = 0;
				y[0] = x[0];

				//Define vertices.
				for(int k=0; k<4; k++)
				{
					oldvertex_y[k] = LATTICE(y[0],y[1],k);
					newvertex_x[k] = oldvertex_x[k];
					newvertex_y[k] = oldvertex_y[k];
				}

				//If no bond exists, create a new one.
				if( oldvertex_x[UP] == NO_BOND )
				{
					newvertex_x[UP] = BOND;
					newvertex_y[DOWN] = BOND;
				}
				//If a bond already exists, delete it.
				else
				{
					newvertex_x[UP] = NO_BOND;
					newvertex_y[DOWN] = NO_BOND;
				}

				break;

			case(LEFT):

				//Define the new location y.
				if( x[0] > 0 )
					y[0] = x[0] - 1;
				else
					y[0] = N0-1;
				y[1] = x[1];

				//Define vertices.
				for(int k=0; k<4; k++)
				{
					oldvertex_y[k] = LATTICE(y[0],y[1],k);
					newvertex_x[k] = oldvertex_x[k];
					newvertex_y[k] = oldvertex_y[k];
				}

				//If no bond exists, create a new one.
				if( oldvertex_x[LEFT] == NO_BOND )
				{
					newvertex_x[LEFT] = BOND;
					newvertex_y[RIGHT] = BOND;
				}
				//If a bond already exists, delete it.
				else
				{
					newvertex_x[LEFT] = NO_BOND;
					newvertex_y[RIGHT] = NO_BOND;
				}

				break;

			case(DOWN):

				//Define the new location y.
				if( x[1] > 0 )
					y[1] = x[1] - 1;
				else
					y[1] = N1-1;
				y[0] = x[0];

				//Define vertices.
				for(int k=0; k<4; k++)
				{
					oldvertex_y[k] = LATTICE(y[0],y[1],k);
					newvertex_x[k] = oldvertex_x[k];
					newvertex_y[k] = oldvertex_y[k];
				}

				//If no bond exists, create a new one.
				if( oldvertex_x[DOWN] == NO_BOND )
				{
					newvertex_x[DOWN] = BOND;
					newvertex_y[UP] = BOND;
				}
				//If a bond already exists, delete it.
				else
				{
					newvertex_x[DOWN] = NO_BOND;
					newvertex_y[UP] = NO_BOND;
				}

				break;

			}
		}

		//Determine the Vertex-Weights.
		double weightx, weighty, weightx_new, weighty_new;

		//Determine the vertex-weights if one removes the worm.
		if(check_end_prop)
		{
			weightx = weight_7Vertex_Cuda(M,oldvertex_x,true);
            weightx_new = weight_7Vertex_Cuda(M,newvertex_x,false);
            weighty = 1.0;
            weighty_new = 1.0;
		}
		else
		{
			//Determine the vertex-weights at x if one is at the start-location.
			if( x[0] == start[0] && x[1] == start[1] )
			{
				weightx = weight_7Vertex_Cuda(M,oldvertex_x,true);

				//Do not accept leaving behind a 3-vertex at the start!
				int n_bonds = 0;
				for(int k=0; k<4; k++)
				{
					if( newvertex_x[k] != 0 )
						++n_bonds;
				}

				if( n_bonds >= 3 )
					weightx_new = 0.0;
				else
					weightx_new = weight_7Vertex_Cuda(M,newvertex_x,true);
			}
			//No special conditions: Worm moves away from x.
			else
			{
				weightx = weight_7Vertex_Cuda(M,oldvertex_x,true);
				weightx_new = weight_7Vertex_Cuda(M,newvertex_x,false);
			}

			//Determine the vertex-weights at y if one moves to the start-location.
			if( y[0] == start[0] && y[1] == start[1] )
			{
				weighty = weight_7Vertex_Cuda(M,oldvertex_y,true);
				weighty_new = weight_7Vertex_Cuda(M,newvertex_y,true);
			}
			//No special conditions: Worm moves to y.
			else
			{
				weighty = weight_7Vertex_Cuda(M,oldvertex_y,false);
				weighty_new = weight_7Vertex_Cuda(M,newvertex_y,true);
			}
		}

		//Determine the proposal and accpetance probability.
		const double q1 = proposal_7Vertex_Cuda(x,y,oldvertex_x,p0,N0,N1);
		const double q2 = proposal_7Vertex_Cuda(y,x,newvertex_y,p0,N0,N1);
		const double pA = (weightx_new/weightx) * (weighty_new/weighty) * (q2/q1);

		//Test if the step will be accepted.
		const double u = randomdouble_Cuda(rng);
		if(u < pA)
		{
			//If one removes the worm.
			if(check_end_prop)
			{
				check_end = true;
				break;
			}
			//If one continues, apply the changes.
			else
			{
				for(int k=0; k<4; k++)
				{
					LATTICE(x[0],x[1],k) = newvertex_x[k];
					LATTICE(y[0],y[1],k) = newvertex_y[k];
				}

				x[0] = y[0];
				x[1] = y[1];
			}
		}
	}
}

__global__ void generateWormData_7Vertex_Cuda(const int* N0_ptr, const int* N1_ptr, const double* M_start_ptr, const int* n_data_ptr, pcg32_random_t* rng)
{
	const int N0 = *N0_ptr, N1 = *N1_ptr, n_data = *n_data_ptr;
	
	//Determine the Mass via the Thread-ID.
	const double M = 0.2 * (double)threadIdx.x + (*M_start_ptr);

	//Initalize Lattice with malloc (CUDA does not support calloc).
	enum BondType* lattice = (enum BondType*)malloc(N0*N1*4*sizeof(int));
	for( int i=0; i<N0*N1*4; i++)
	{
		lattice[i] = NO_BOND;
	}

	//Generate n_data-1 more closed loops.
	for(int i=1; i<n_data;i++)
	{
		Worm_7Vertex_Cuda(N0,N1,M,lattice,rng);
	}

	free(lattice);
}

void cuda7VertexWorm(int N0, int N1, const double M_start, const int n_data, const int n_mass)
{
	//GPU-Variables.
	int *N0d, *N1d, *n_data_d;
	double *M_start_d;
	pcg32_random_t *rng_d;

	//Allocate space on the GPU.
	cudaMalloc((void**)&N0d,sizeof(int));
	cudaMalloc((void**)&N1d,sizeof(int));
	cudaMalloc((void**)&n_data_d,sizeof(int));
	cudaMalloc((void**)&M_start_d,sizeof(double));
	cudaMalloc((void**)&rng_d,sizeof(pcg32_random_t));

	//Copy data to the GPU.
	cudaMemcpy( N0d, &N0, sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy( N1d, &N1, sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy( n_data_d, &n_data, sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy( M_start_d, &M_start, sizeof(double), cudaMemcpyHostToDevice );
	cudaMemcpy( rng_d, &pcg32_global, sizeof(int), cudaMemcpyHostToDevice );

	//Start Kernels.
	generateWormData_7Vertex_Cuda<<<1,n_mass>>>(N0d,N1d,M_start_d,n_data_d,rng_d);

	//Wait for the Kernels to finish.
	cudaDeviceSynchronize();

	//Free Memory.
	cudaFree(N0d);
	cudaFree(N1d);
	cudaFree(M_start_d);
	cudaFree(n_data_d);
	cudaFree(rng_d);
}