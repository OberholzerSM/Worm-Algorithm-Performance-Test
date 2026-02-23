#include <iostream>
#include <vector>
#include <thread>
#include "FortranInterop.h"
#include "TimerHeader.h"
#include "Worm7Vertex_C.h"
#include "Worm7Vertex_CPP.h"
#include "RandomHeader.h"
#include "Worm7Vertex_CUDA.cuh"

Timer CLOCK{};

enum class Option
{
	TEST_FORTRAN_VERSIONS,
	TEST_FORTRAN_VERSIONS_SHORT,
	TEST_LANGUAGE_VERSIONS,
	TEST_RNG_VERSIONS_CPP,
	TEST_PARALLEL_VERSIONS,
};

void testFortranVersions(int N0, int N1, const float M_start = 0.0f, const int n_data = (int)10E5, const int n_mass = 21)
{
	constexpr int n_versions = 10;
	float M;
	double t1, t2;
	void (*p[n_versions]) (const int *N0, const int *N1, const float *M, const int *n_data) = 
	{ generate7VertexWormDataFortran_v1, generate7VertexWormDataFortran_v2, generate7VertexWormDataFortran_v3, generate7VertexWormDataFortran_v4, generate7VertexWormDataFortran_v5,
	  generate7VertexWormDataFortran_v6, generate7VertexWormDataFortran_v7, generate7VertexWormDataFortran_v8, generate7VertexWormDataFortran_v9, generate7VertexWormDataFortran_v10};

	std::cout << "N=" << N0 << 'x' << N1 << '\n';
	for(int n=0; n<n_versions; n++)
	{
		M = M_start;
		t1 = CLOCK.getTime();
		for(int i=0; i<n_mass; i++)
		{
			p[n](&N0, &N1, &M, &n_data);
			M += 0.2f;
		}
		t2 = CLOCK.getTime();
		std::cout << "Time Fortran v" << n+1 << ": " << t2-t1 << "s\n";
	}
	std::cout << '\n';
}

void testFortranVersions_Short(int N0, int N1, const float M_start = 0.0f, const int n_data = (int)10E5, const int n_mass = 21)
{
	constexpr int n_versions = 4;
	float M;
	double t1, t2;
	void (*p[n_versions]) (const int *N0, const int *N1, const float *M, const int *n_data) =
	{generate7VertexWormDataFortran_v2, generate7VertexWormDataFortran_v3, generate7VertexWormDataFortran_v9, generate7VertexWormDataFortran_v10};
	int version[n_versions] = {2,3,9,10};

	std::cout << "N=" << N0 << 'x' << N1 << '\n';
	for(int n=0; n<n_versions; n++)
	{
		M = M_start;
		t1 = CLOCK.getTime();
		for(int i=0; i<n_mass; i++)
		{
			p[n](&N0, &N1, &M, &n_data);
			M += 0.2f;
		}
		t2 = CLOCK.getTime();
		std::cout << "Time Fortran v" << version[n] << ": " << t2-t1 << "s\n";
	}
	std::cout << '\n';
}

void testLanguageVersions(int N0, int N1, const double M_start = 0.0, const int n_data = (int)10E5, const int n_mass = 21)
{
	double M;
	float Mf;
	double t1, t2;
	std::cout << "N=" << N0 << 'x' << N1 << '\n';

	//Fortran Version with the best time.
	Mf = (float)M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		generate7VertexWormDataFortran_v10(&N0,&N1,&Mf,&n_data);
		Mf += 0.2f;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time Fortran:\t" << t2-t1 << "s\n";

	//C Version
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		generateWormData_7Vertex_C(N0,N1,M,n_data);
		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time C:\t\t" << t2-t1 << "s\n";

	//C++ Version
	Random::mode = Random::PCG32;
	Random::setSeed();
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time C++:\t" << t2-t1 << "s\n";

	std::cout << '\n';
}

void testRNGVersionsCPP(int N0, int N1, const double M_start = 0.0, const int n_data = (int)10E5, const int n_mass = 21)
{
	double M;
	double t1, t2;
	std::cout << "N=" << N0 << 'x' << N1 << '\n';

	//C++ Version using C RNG
	Random::mode = Random::RAND;
	Random::setSeed();
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time RAND:\t" << t2-t1 << "s\n";

	Random::mode = Random::MTWISTER32;
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time Twister32:\t" << t2-t1 << "s\n";

	Random::mode = Random::MTWISTER64;
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time Twister64:\t" << t2-t1 << "s\n";

	Random::mode = Random::RANLUX24;
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time Ranlux24:\t" << t2-t1 << "s\n";

	Random::mode = Random::RANLUX48;
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time Ranlux48:\t" << t2-t1 << "s\n";

	Random::mode = Random::SPLITMIX64;
	Random::setSeed();
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time SplitMix:\t" << t2-t1 << "s\n";

	Random::mode = Random::PCG32;
	Random::setSeed();
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time PCG:\t" << t2-t1 << "s\n";

	Random::mode = Random::XORSHIFTM;
	Random::setSeed();
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time xorshift*:\t" << t2-t1 << "s\n";

	Random::mode = Random::XOSHIRO256MM;
	Random::setSeed();
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		//Create Worm object.
		Worm<DiracSCNf1> worm{N0,N1,M};

		//Do n_data many loops.
		for(int j=1; j<n_data; j++)
		{
			worm.doClosedLoop();
		}

		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time xoshiro256**:\t" << t2-t1 << "s\n";

	std::cout << '\n';
}

void testParalellVersions(int N0, int N1, const double M_start = 0.0, const int n_data = (int)10E5, const int n_mass = 21)
{
	double M;
	double t1, t2;
	std::cout << "N=" << N0 << 'x' << N1 << '\n';

	//Default C version.
	M = M_start;
	t1 = CLOCK.getTime();
	for(int i=0; i<n_mass; i++)
	{
		generateWormData_7Vertex_C(N0,N1,M,n_data);
		M += 0.2;
	}
	t2 = CLOCK.getTime();
	std::cout << "Time Default:\t" << t2-t1 << "s\n";

	//Thread-Version.
	t1 = CLOCK.getTime();
	std::vector<std::thread> threadList{};

	//Launch Threads
	M = M_start;
	for(int i=0; i<n_mass; i++)
	{
		threadList.emplace_back( [N0,N1,M,n_data]() {generateWormData_7Vertex_C(N0,N1,M,n_data);} );
		M += 0.2;
	}

	//Wait for the Threads to finish.
	for(int i=0; i<n_mass; i++)
	{
		threadList[i].join();
	}

	t2 = CLOCK.getTime();
	std::cout << "Time Threads:\t" << t2-t1 << "s\n";

	//CUDA-Version (slow, only do if the volume is small)
	if(N0*N1<=4)
	{
		t1 = CLOCK.getTime();
		cuda7VertexWorm(N0,N1,M_start,n_data,n_mass);
		t2 = CLOCK.getTime();
		std::cout << "Time CUDA:\t" << t2-t1 << "s\n";
	}
	
	std::cout << '\n';
}

int main(int argc, char* argv)
{
	constexpr Option option = Option::TEST_PARALLEL_VERSIONS;
	switch(option)
	{
	case Option::TEST_FORTRAN_VERSIONS:
		testFortranVersions(2, 2);
		break;
		//Result: No significant difference between v2-v9 for 2x2 lattices. v10 slightly faster.

	case Option::TEST_FORTRAN_VERSIONS_SHORT:
		testFortranVersions_Short(2, 2);
		testFortranVersions_Short(32, 2);
		break;
		//Result: v10 is the fastest, with v2 being second place.

	case Option::TEST_LANGUAGE_VERSIONS:
		testLanguageVersions(2,2);
		testLanguageVersions(32,2);
		break;
		//Result: Using PCG, C appears to be the fastest and Fortran the slowest.

	case Option::TEST_RNG_VERSIONS_CPP:
		testRNGVersionsCPP(2,2);
		break;
		//Result: PCG and xorshift* are the fastest.

	case Option::TEST_PARALLEL_VERSIONS:
		testParalellVersions(2,2);
		testParalellVersions(32,2);
		//Result: Threads make it faster on small latices but slower on larger ones, CUDA makes it significantly slower.
		break;
	}

	return 0;
}