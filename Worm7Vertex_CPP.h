#pragma once
#include <iostream>
#include <string>
#include <stdlib.h>

//Represents a Point on the Lattice in Lattice-Coordinates.
struct Point2D
{
	int i,j;

	Point2D& operator=(const Point2D &point);
	Point2D& operator=(const int x);
};
std::ostream &operator<<(std::ostream &out, const Point2D &point);
bool operator== (const Point2D& p1, const Point2D& p2);

//Possible Directions on the Lattice
enum Direction
{
	RIGHT,
	UP,
	LEFT,
	DOWN,
};

std::string directionToString(Direction direction);

//Worm-Types

enum class WormType
{
	DIRACSCNF1,
};

//Bond-Types for each Worm-Type

//Strong Coupling Dirac with Nf=1.
enum class DiracSCNf1
{
	NONE,
	BOND,
};

//Vertex-Struct.
template <typename Bond>
struct Vertex
{
	Bond mu[4] = {Bond::NONE,Bond::NONE,Bond::NONE,Bond::NONE};
	int getBonds() const;
};

//Worm on a Lattice.
template <typename Bond>
class Worm
{
public:

	Worm(int N0 = 2, int N1 = 2, double M = 0.0) : lattice{lattice}, N0{N0}, N1{N1}, M{M}
	{
		lattice = new Vertex<Bond>[N0 * N1];
		//p0 = 1.0 / (double)(N0*N1);
		p0 = 0.5;
		srand((unsigned int)time(NULL)); 
	}

	//Create a new closed loop.
	void doClosedLoop();

	//Return the Vertex-Weight.
	double vertexWeight(const Vertex<Bond> &vertex, bool sink = false) const;

	//Returns the probability of proposing a change.
	double proposalProb(const Point2D &loc_old, const Point2D &loc_new, const Vertex<Bond> &oldvertex_x) const;

	//Returns the vertex at a given point.
	Vertex<Bond> getVertex(const Point2D &point) const;

	~Worm()
	{
		delete[] lattice;
	}

	bool debug = false;

private:

	int N0,N1;							//Lattice Size.
	double M;							//Bare Mass parameter.
	double p0{};						//Probability to remove the Worm from the lattice.
	Vertex<Bond>* lattice;				//Lattice of Vertices.
	Point2D x{}, y{}, start{};			//Location of the wormhead, proposed new location for the wormhead and the starting location.
};