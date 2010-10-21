/*
 * PolyLine_Simplification.h
 *
 *  Created on: Aug 15, 2010
 *      Author: ptroja
 */

#ifndef POLYLINE_SIMPLIFICATION_H_
#define POLYLINE_SIMPLIFICATION_H_

#include <vector>
#include <valarray>

//===================================================================
// Definitions
//===================================================================
typedef double Scalar;

typedef std::valarray<double> PointBase;

class Point : public PointBase
{
public:
	Point();
	Point(const PointBase & va);
	Point(double _x, double _y);

    // Reuse assignment operators from base class
    using PointBase::operator=;
};

class Vector : public PointBase
{
public:
	Vector();
	Vector(const PointBase & va);

	// Reuse assignment operators from base class
	using PointBase::operator=;
};

class Segment
{
public:
	Point P0, P1;
};

// dot product (3D) which allows vector operations in arguments
#define dot(u,v)   (((u)*(v)).sum())
#define norm2(v)   dot((v),(v))        // norm2 = squared length of vector
#define norm(v)    std::sqrt(norm2(v))      // norm = length of vector
#define d2(u,v)    norm2((u)-(v))      // distance squared = norm2 of difference
#define d(u,v)     norm((u)-(v))       // distance = norm of difference

void simplifyDP(const double tol, const Point* v, const int j, const int k, bool * mk);
unsigned int poly_simplify(const double tol, const std::vector<Point> & V, std::vector<Point> & sV);
double dist_Point_to_Segment(const Point & P, const Segment & S, Point & C);

#endif /* POLYLINE_SIMPLIFICATION_H_ */
