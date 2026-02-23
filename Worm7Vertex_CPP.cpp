#include "Worm7Vertex_CPP.h"
#include "RandomHeader.h"
#include <iostream>

//Operator Overloading
std::ostream &operator<<(std::ostream &out, const Point2D &point)
{
	out << point.i << ',' << point.j;
	return out;
}

Point2D& Point2D::operator=(const Point2D &point)
{
	i = point.i;
	j = point.j;
	return *this;
}

Point2D& Point2D::operator=(const int x)
{
	i = x;
	j = x;
	return *this;
}

bool operator== (const Point2D& p1, const Point2D& p2)
{
	return ( (p1.i == p2.i) && (p1.j == p2.j) );
}

//Enum to String
std::string directionToString(Direction direction)
{
	switch(direction)
	{
	case RIGHT:
		return "right";

	case UP:
		return "up";

	case LEFT:
		return "left";

	case DOWN:
		return "down";

	default:
		return "?";
	}
}

//Vertex-Function.
template <typename Bond>
int Vertex<Bond>::getBonds() const
{
	int nBonds = 0;
	for(int k=0; k<4; k++)
	{
		if( mu[k] != Bond::NONE )
			++nBonds;
	}
	return nBonds;
}

//Worm Functions

template <typename Bond>
Vertex<Bond> Worm<Bond>::getVertex(const Point2D &point) const
{
	return lattice[ (point.i-1)*N1 + (point.j-1) ];
}

//Dirac Strong Coupling Nf=1 Worm functions.

void Worm<DiracSCNf1>::doClosedLoop()
{
	using enum DiracSCNf1;
	using namespace Random;

	//Propose a new start.
	start.i = randomint(1,N0);
	start.j = randomint(1,N1);

	//Determine the acceptance probability.
	const Vertex<DiracSCNf1> vertexStart = getVertex(start);
	const double weight_old = vertexWeight(vertexStart,false);
	const double weight_new = vertexWeight(vertexStart,true);
	const double q1 = 1.0 / (double)(N0*N1);
	const double q2 = p0;
	const double pA = (weight_new / weight_old) * (q2/q1);

	//Test if the new start will be accepted.
	const double u = randomdouble();
	if(u < pA)
		x = start;
	else
		return;

	while(1)
	{
		//Set the current vertex.
		const Vertex<DiracSCNf1> oldvertex_x = getVertex(x);

		//Declare other vertices.
		Vertex<DiracSCNf1> newvertex_x, oldvertex_y, newvertex_y;

		//If a loop was completed, propose to remove the worm.
		bool check_end_prop = false;
		if(x == start)
		{
			const double u = randomdouble();
			if(u < p0)
			{
				check_end_prop = true;
				y = 0;
				for(int k=0; k<4; k++)
				{
					newvertex_x.mu[k] = oldvertex_x.mu[k];
					oldvertex_y.mu[k] = oldvertex_x.mu[k];
					newvertex_y.mu[k] = oldvertex_x.mu[k];
				}
			}
		}

		//If the worm will not be removed, propose a new location.
		if(!check_end_prop)
		{
			//Propose a new direction to move to.
			const Direction direction = (Direction)randomint(0,3);

			switch(direction)
			{
			case RIGHT:
				//Define the new location y.
				if( x.i < N0 )
					y.i = x.i + 1;
				else
					y.i = 1;
				y.j = x.j;
				oldvertex_y = getVertex(y);

				//Define Vertices.
				for(int k=0; k<4; k++)
				{
					newvertex_x.mu[k] = oldvertex_x.mu[k];
					newvertex_y.mu[k] = oldvertex_y.mu[k];
				}

				//If no bond exists, add one.
				if( oldvertex_x.mu[direction] == NONE )
				{
					newvertex_x.mu[RIGHT] = BOND;
					newvertex_y.mu[LEFT] = BOND;
				}
				else
				{
					newvertex_x.mu[RIGHT] = NONE;
					newvertex_y.mu[LEFT] = NONE;
				}

				break;

			case UP:
				//Define the new location y.
				if( x.j < N1 )
					y.j = x.j + 1;
				else
					y.j = 1;
				y.i = x.i;
				oldvertex_y = getVertex(y);

				//Define Vertices.
				for(int k=0; k<4; k++)
				{
					newvertex_x.mu[k] = oldvertex_x.mu[k];
					newvertex_y.mu[k] = oldvertex_y.mu[k];
				}

				//If no bond exists, add one.
				if( oldvertex_x.mu[direction] == NONE )
				{
					newvertex_x.mu[UP] = BOND;
					newvertex_y.mu[DOWN] = BOND;
				}
				else
				{
					newvertex_x.mu[UP] = NONE;
					newvertex_y.mu[DOWN] = NONE;
				}

				break;

			case LEFT:
				//Define the new location y.
				if( x.i > 1 )
					y.i = x.i - 1;
				else
					y.i = N0;
				y.j = x.j;
				oldvertex_y = getVertex(y);

				//Define Vertices.
				for(int k=0; k<4; k++)
				{
					newvertex_x.mu[k] = oldvertex_x.mu[k];
					newvertex_y.mu[k] = oldvertex_y.mu[k];
				}

				//If no bond exists, add one.
				if( oldvertex_x.mu[direction] == NONE )
				{
					newvertex_x.mu[LEFT] = BOND;
					newvertex_y.mu[RIGHT] = BOND;
				}
				else
				{
					newvertex_x.mu[LEFT] = NONE;
					newvertex_y.mu[RIGHT] = NONE;
				}

				break;

			case DOWN:
				//Define the new location y.
				if( x.j > 1 )
					y.j = x.j - 1;
				else
					y.j = N1;
				y.i = x.i;
				oldvertex_y = getVertex(y);

				//Define Vertices.
				for(int k=0; k<4; k++)
				{
					newvertex_x.mu[k] = oldvertex_x.mu[k];
					newvertex_y.mu[k] = oldvertex_y.mu[k];
				}

				//If no bond exists, add one.
				if( oldvertex_x.mu[direction] == NONE )
				{
					newvertex_x.mu[DOWN] = BOND;
					newvertex_y.mu[UP] = BOND;
				}
				else
				{
					newvertex_x.mu[DOWN] = NONE;
					newvertex_y.mu[UP] = NONE;
				}

				break;
			}
		}

		//Determine the Vertex-Weights.
		double weightx, weighty, weightx_new, weighty_new;

		if(check_end_prop)
		{
			weightx = vertexWeight(oldvertex_x,true);
			weightx_new = vertexWeight(newvertex_x,false);
			weighty = 1.0;
			weighty_new = 1.0;
		}
		else
		{
			//Determine the Vertex-Weights at x.
			if( x == start )
			{
				weightx = vertexWeight(oldvertex_x,true);
				const int nBonds = newvertex_x.getBonds();
				weightx_new = ( nBonds == 3 ? 0.0 : vertexWeight(newvertex_x,true) );
			}
			else
			{
				weightx = vertexWeight(oldvertex_x,true);
				weightx_new = vertexWeight(newvertex_x,false);
			}

			//Determine the Vertex-Weights at y.
			if( y == start )
			{
				weighty = vertexWeight(oldvertex_y,true);
				weighty_new = vertexWeight(newvertex_y,true);
			}
			else
			{
				weighty = vertexWeight(oldvertex_y,false);
				weighty_new = vertexWeight(newvertex_y,true);
			}
		}

		//Determine proposal and acceptance probabilities.
		const double q1 = proposalProb(x,y,oldvertex_x);
		const double q2 = proposalProb(y,x,newvertex_y);
		const double pA = (weightx_new/weightx) * (weighty_new/weighty) * (q2/q1);

		//Test if the step will be accepted.
		const double u = randomdouble();
		if(u < pA)
		{
			//If one removes the worm.
			if(check_end_prop)
				return;

			for(int k=0; k<4; k++)
			{
				lattice[ (x.i-1)*N1 + (x.j-1) ].mu[k] = newvertex_x.mu[k];
				lattice[ (y.i-1)*N1 + (y.j-1) ].mu[k] = newvertex_y.mu[k];
			}
			x = y;
		}
	}
}

template <>
double Worm<DiracSCNf1>::vertexWeight(const Vertex<DiracSCNf1> &vertex, bool sink) const
{
	constexpr double ROOT3_INV4 = 0.6299605249474366;
	const int nBonds = vertex.getBonds();	
	using enum DiracSCNf1;
	const bool straightLine = ( vertex.mu[RIGHT] == BOND && vertex.mu[LEFT] == BOND ) || ( vertex.mu[UP] == BOND && vertex.mu[DOWN] == BOND );

	switch(nBonds)
	{
	case 0:
		return sink ? 4.0 : M*M;

	case 1:
		return sink ? 1.0 : 0.0;

	case 2:
		return straightLine ? 1.0 : 0.5;

	case 3:
		return sink ? ROOT3_INV4 : 0.0;

	default:
		return 0.0;
	}
}

template <>
double Worm<DiracSCNf1>::proposalProb(const Point2D &loc_old, const Point2D &loc_new, const Vertex<DiracSCNf1> &oldvertex_x) const
{
	//If outside the Worm.
	if( loc_old.i == 0 )
		return 1.0 / (double)(N0*N1);

	const int nBonds = oldvertex_x.getBonds();

	switch(nBonds)
	{
	case 0:
	case 2:
		return ( loc_new.i == 0 ) ? p0 : (1.0-p0)/4.0;

	case 1:
	case 3:
		return 0.25;

	default:
		return 0.0;
	}
}