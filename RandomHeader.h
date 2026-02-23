#pragma once
#include <stdlib.h>
#include <time.h>
#include <random> 

namespace Random
{
	//What RNG to use.
	enum RNG
	{
		RAND,
		MTWISTER32,
		MTWISTER64,
		RANLUX24,
		RANLUX48,
		PCG32,
		SPLITMIX64,
		XORSHIFTM,
		XOSHIRO256MM,
	};
	inline RNG mode = RANLUX24;

	//C++
	inline auto u_distribution = std::uniform_real_distribution{ 0.0, 1.0 };			//To create a random double between 0 and 1.

	//Mersenne-Twister.
	namespace MT32
	{
		inline std::random_device rd{};													//random_device: Random seed determined by the computer.
		inline std::seed_seq seed{ rd(),rd(), rd(), rd(), rd(), rd(), rd(), rd() };		//Generate a new seed using eight seeds from random_device.
		inline std::mt19937 generator{ seed };											//Globale mt19937 object. Generates a 32it integer using the seed.
	}

	namespace MT64
	{
		inline std::random_device rd{};													//random_device: Random seed determined by the computer.
		inline std::seed_seq seed{ rd(),rd(), rd(), rd(), rd(), rd(), rd(), rd() };		//Generate a new seed using eight seeds from random_device.
		inline std::mt19937_64 generator{ seed };										//Globale mt19937 object. Generates a 64it integer using the seed.
	}

	//Ranlux
	namespace Ranlux24
	{
		inline std::random_device rd{};													//random_device: Random seed determined by the computer.
		inline std::seed_seq seed{ rd(),rd(), rd(), rd(), rd(), rd(), rd(), rd() };		//Generate a new seed using eight seeds from random_device.
		inline std::ranlux24 generator{ seed };		
	}

	namespace Ranlux48
	{
		inline std::random_device rd{};													//random_device: Random seed determined by the computer.
		inline std::seed_seq seed{ rd(),rd(), rd(), rd(), rd(), rd(), rd(), rd() };		//Generate a new seed using eight seeds from random_device.
		inline std::ranlux48 generator{ seed };		
	}

	//SplitMix64
	namespace splitmix64
	{
		inline uint64_t x = 1;

		inline void setSeed(uint64_t seed)
		{
			x = seed;
		}

		inline uint64_t get64uint()
		{
			uint64_t z = (x += 0x9e3779b97f4a7c15);
			z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
			z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
			return z ^ (z >> 31);
		}
	}

	//Permuted Congruential Generator
	namespace pcg32
	{
		struct pcg_state_setseq_64 {    // Internals are *Private*.
			uint64_t state;             // RNG state.  All values are possible.
			uint64_t inc;               // Controls which RNG sequence (stream) is
										// selected. Must *always* be odd.
		};
		typedef struct pcg_state_setseq_64 pcg32_random_t;

		// If you *must* statically initialize it, here's one.

		#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }
		static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

		// pcg32_random()
		// pcg32_random_r(rng)
		//     Generate a uniformly distributed 32-bit random number

		inline uint32_t pcg32_random_r(pcg32_random_t* rng)
		{
			uint64_t oldstate = rng->state;
			rng->state = oldstate * 6364136223846793005ULL + rng->inc;
			uint32_t xorshifted = (uint32_t) ( ((oldstate >> 18u) ^ oldstate) >> 27u );
			uint32_t rot = (uint32_t)(oldstate >> 59u);
			return (xorshifted >> rot) | (xorshifted << ((0-rot) & 31));
		}

		inline uint32_t pcg32_random(void)
		{
			 return pcg32_random_r(&pcg32_global);
		}

		// pcg32_srandom(initstate, initseq)
		// pcg32_srandom_r(rng, initstate, initseq):
		//     Seed the rng.  Specified in two parts, state initializer and a
		//     sequence selection constant (a.k.a. stream id)

		inline void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
		{
			rng->state = 0U;
			rng->inc = (initseq << 1u) | 1u;
			pcg32_random_r(rng);
			rng->state += initstate;
			pcg32_random_r(rng);
		}

		inline void pcg32_srandom(uint64_t seed, uint64_t seq)
		{
			pcg32_srandom_r(&pcg32_global, seed, seq);
		}


		// pcg32_boundedrand(bound):
		// pcg32_boundedrand_r(rng, bound):
		//     Generate a uniformly distributed number, r, where 0 <= r < bound

		inline uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)
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

		inline uint32_t pcg32_boundedrand(uint32_t bound)
		{
			return pcg32_boundedrand_r(&pcg32_global, bound);
		}
	}

	//xorshift*
	namespace xorshiftp
	{
		inline uint64_t x = 1;

		inline void setSeed(uint64_t seed)
		{
			x = seed;
		}

		inline uint32_t get32uint(void)
		{
			x ^= x >> 12;
			x ^= x << 25;
			x ^= x >> 27;
			return (uint32_t)(x * 0x2545F4914F6CDD1DULL);
		}
	}

	//xoshiro256**
	namespace xoshiro256pp
	{
		//Helper-Function.
		inline uint64_t rol64(uint64_t x, int k) 
		{
			return (x << k) | (x >> (64 - k));
		}

		inline uint64_t s[4];

		inline void setSeed(uint64_t seed)
		{
			splitmix64::setSeed(seed);
			for(int i=0; i<4; i++)
			{
				s[i] = splitmix64::get64uint();
			}
		}

		inline uint64_t get64uint() 
		{
			const uint64_t result = rol64(s[1] * 5, 7) * 9;
			const uint64_t t = s[1] << 17;

			s[2] ^= s[0];
			s[3] ^= s[1];
			s[1] ^= s[2];
			s[0] ^= s[3];

			s[2] ^= t;
			s[3] = rol64(s[3], 45);

			return result;
		}
	}

	//Set the Seed.
	inline void setSeed()
	{
		unsigned int t = (unsigned int)time(NULL);
		switch(mode)
		{
		case RAND:
			srand(t);
			return;

		//C++ PRNGS: Seed already set via object creation.
		case MTWISTER32:
		case MTWISTER64:
		case RANLUX24:
		case RANLUX48:
			return;

		case SPLITMIX64:
			splitmix64::setSeed(t);
			return;

		case PCG32:
			//pcg32::setSeed(t);
			return;

		case XORSHIFTM:
			xorshiftp::setSeed(t);
			return;

		case XOSHIRO256MM:
			xoshiro256pp::setSeed(t);
			return;

		default:
			return;
		}
	}

	//Return a random int based on the chosen PRNG.
	inline int randomint(const int a, const int b)
	{
		switch(mode)
		{
		case RAND:
		{
			const int d = b-a+1;
			const int end = (RAND_MAX/d)*d;
			int r = rand();
			while(r >= end)
				r = rand();	
			return r%d + a;
		}

		case MTWISTER32:
			return std::uniform_int_distribution{ a, b }(MT32::generator);

		case MTWISTER64:
			return std::uniform_int_distribution{ a, b }(MT64::generator);

		case RANLUX24:
			return std::uniform_int_distribution{ a, b }(Ranlux24::generator);

		case RANLUX48:
			return std::uniform_int_distribution{ a, b }(Ranlux48::generator);

		case SPLITMIX64:
		{
			const uint64_t d = (uint64_t)(b-a+1);
			const uint64_t end = (UINT64_MAX/d)*d;
			uint64_t r = splitmix64::get64uint();
			while(r >= end)
				r = splitmix64::get64uint();
			return (int)(r%d + a);
		}

		case PCG32:
		{
			const uint32_t d = (uint32_t)(b-a+1);
			const uint32_t r = pcg32::pcg32_boundedrand(d);
			return (int)r + a;
		}

		case XORSHIFTM:
		{
			const uint32_t d = (uint32_t)(b-a+1);
			const uint32_t end = (UINT32_MAX/d)*d;
			uint32_t r = xorshiftp::get32uint();
			while(r >= end)
				r = xorshiftp::get32uint();
			return (int)(r%d + a);
		}

		case XOSHIRO256MM:
		{
			const uint64_t d = (uint64_t)(b-a+1);
			const uint64_t end = (UINT64_MAX/d)*d;
			uint64_t r = xoshiro256pp::get64uint();
			while(r >= end)
				r = xoshiro256pp::get64uint();
			return (int)(r%d + a);
		}

		default:
			return std::uniform_int_distribution{ a, b }(MT64::generator);

		}
	}

	//Return a random double based on the chosen PRNG.
	inline double randomdouble(const double min = 0.0, const double max = 1.0)
	{
		switch(mode)
		{
		case RAND:
		{
			const double u = (double)rand() / (double)RAND_MAX;
			return u*(max-min) + min;
		}

		case MTWISTER32:
			if( min == 0.0 && max == 1.0 )
				return u_distribution(MT32::generator);
			else
				return std::uniform_real_distribution{ min, max }(MT32::generator);

		case MTWISTER64:
			if( min == 0.0 && max == 1.0 )
				return u_distribution(MT64::generator);
			else
				return std::uniform_real_distribution{ min, max }(MT64::generator);

		case RANLUX24:
			if( min == 0.0 && max == 1.0 )
				return u_distribution(Ranlux24::generator);
			else
				return std::uniform_real_distribution{ min, max }(Ranlux24::generator);

		case RANLUX48:
			if( min == 0.0 && max == 1.0 )
				return u_distribution(Ranlux48::generator);
			else
				return std::uniform_real_distribution{ min, max }(Ranlux48::generator);

		case SPLITMIX64:
		{
			const double u = (double)splitmix64::get64uint() / (double)UINT64_MAX;
			return u*(max-min) + min;
		}

		case PCG32:
		{
			const double u = (double)pcg32::pcg32_random() / (double)UINT32_MAX;
			return u*(max-min) + min;
		}

		case XORSHIFTM:
		{
			const double u = (double)xorshiftp::get32uint() / (double)UINT32_MAX;
			return u*(max-min) + min;
		}

		case XOSHIRO256MM:
		{
			const double u = (double)xoshiro256pp::get64uint() / (double)UINT64_MAX;
			return u*(max-min) + min;
		}
		
		default:
			if( min == 0.0 && max == 1.0 )
				return u_distribution(MT64::generator);
			else
				return std::uniform_real_distribution{ min, max }(MT64::generator);

		}
	}
}