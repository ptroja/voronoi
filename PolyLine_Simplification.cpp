// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// Assume that classes are already given for the objects:
//    Point and Vector with
//        coordinates {double x, y;}    // as many as are needed
//        operators for:
//            == to test equality
//            != to test inequality
//            (Vector)0 = (0,0,0)         (null vector)
//            Point  = Point ± Vector
//            Vector = Point - Point
//            Vector = Vector ± Vector
//            Vector = Scalar * Vector    (scalar product)
//            Vector = Vector * Vector    (cross product)
//    Segment with defining endpoints {Point P0, P1;}
//===================================================================

#include <cmath>
#include <vector>

#include "PolyLine_Simplification.h"

//===================================================================
// Implementation
//===================================================================
Point::Point()
	: PointBase (0.0f, 2)
{
}

Point::Point(const PointBase & va)
	: PointBase (va)
{
}

Point::Point(double _x, double _y)
	: PointBase (2)
{
	this->operator[](0) = _x;
	this->operator[](1) = _y;
}

Vector::Vector()
	: PointBase (0.0f, 2)
{
}

Vector::Vector(const PointBase & va)
	: PointBase (va)
{
}

// simplifyDP():
//  This is the Douglas-Peucker recursive simplification routine
//  It just marks vertices that are part of the simplified polyline
//  for approximating the polyline subchain v[j] to v[k].
//    Input:  tol = approximation tolerance
//            v[] = polyline array of vertex points
//            j,k = indices for the subchain v[j] to v[k]
//    Output: mk[] = array of markers matching vertex array v[]
void simplifyDP(const double tol, const Point* v, const int j, const int k, bool * mk)
{
	if (k <= j + 1) // there is nothing to simplify
		return;

	// check for adequate approximation by segment S from v[j] to v[k]
	int maxi = j; // index of vertex farthest from S
	double maxd2 = 0; // distance squared of farthest vertex
	const double tol2 = tol * tol; // tolerance squared
	const Segment S = { v[j], v[k] }; // segment from v[j] to v[k]
	const Vector u = Vector(S.P1 - S.P0); // segment direction vector
	const double cu = dot(u,u); // segment length squared

	// test each vertex v[i] for max distance from S
	// compute using the Feb 2001 Algorithm's dist_Point_to_Segment()
	// Note: this works in any dimension (2D, 3D, ...)
	Vector w;
	Point Pb; // base of perpendicular from v[i] to S
	double b, cw, dv2; // dv2 = distance v[i] to S squared

	for (int i = j + 1; i < k; i++) {
		// compute distance squared
		w = v[i] - S.P0;
		cw = dot(w,u);
		if (cw <= 0)
			dv2 = d2(v[i], S.P0);
		else if (cu <= cw)
			dv2 = d2(v[i], S.P1);
		else {
			b = cw / cu;
			Pb = S.P0 + u * b;
			dv2 = d2(v[i], Pb);
		}
		// test with current max distance squared
		if (dv2 <= maxd2)
			continue;
		// v[i] is a new max vertex
		maxi = i;
		maxd2 = dv2;
	}
	if (maxd2 > tol2) // error is worse than the tolerance
	{
		// split the polyline at the farthest vertex from S
		mk[maxi] = true; // mark v[maxi] for the simplified polyline
		// recursively simplify the two subpolylines at v[maxi]
		simplifyDP(tol, v, j, maxi, mk); // polyline v[j] to v[maxi]
		simplifyDP(tol, v, maxi, k, mk); // polyline v[maxi] to v[k]
	}
	// else the approximation is OK, so ignore intermediate vertices
	return;
}

// poly_simplify():
//    Input:  tol = approximation tolerance
//            V[] = polyline array of vertex points
//    Output: sV[]= simplified polyline vertices (max is n)
//    Return: m   = the number of points in sV[]
unsigned int poly_simplify(const double tol, const std::vector<Point> & V, std::vector<Point> & sV)
{
	const unsigned int n = V.size();
	unsigned int i, k, pv; // misc counters
	const double tol2 = tol * tol; // tolerance squared
	Point vt[n]; // vertex buffer
	bool mk[n]; // = {0};  // marker buffer

	// STAGE 1.  Vertex Reduction within tolerance of prior vertex cluster
	vt[0] = V[0]; // start at the beginning
	for (i = k = 1, pv = 0; i < n; i++) {
		if (d2(V[i], V[pv]) < tol2)
			continue;
		vt[k++] = V[i];
		pv = i;
	}
	if (pv < n - 1)
		vt[k++] = V[n - 1]; // finish at the end

	// STAGE 2.  Douglas-Peucker polyline simplification
	mk[0] = mk[k - 1] = true; // mark the first and last vertices
	for(unsigned int j = 1; j < k - 1; ++j) {
		mk[j] = false;
	}
	simplifyDP(tol, vt, 0, k - 1, mk);

	// copy marked vertices to the output simplified polyline
	for (i = 0; i < k; i++) {
		if (mk[i]) {
			sV.push_back(vt[i]);
		}
	}

	return sV.size(); // m vertices in simplified polyline
}

// dist_Point_to_Segment(): get the distance of a point to a segment.
//    Input:  a Point P and a Segment S (in any dimension)
//    Output: a Point C on the S, nearest to P
//    Return: the shortest distance from P to S
double dist_Point_to_Segment(const Point & P, const Segment & S, Point & C)
{
	const Vector v = Vector(S.P1 - S.P0);
	const Vector w = Vector(P - S.P0);

	const double c1 = dot(w,v);
	if (c1 <= 0) {
		C = S.P0;
		return d(P, S.P0);
	}

	const double c2 = dot(v,v);
	if (c2 <= c1) {
		C = S.P1;
		return d(P, S.P1);
	}

	const double b = c1 / c2;
	C = Point(S.P0 + b * v);

	return d(P, C);
}
//===================================================================
