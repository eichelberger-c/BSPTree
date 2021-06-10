#pragma once
#include"stdafx.h"
#include "Object.h"
#include "Node.h"
#include "Polygon.h"

#define MIN_POLY_THRESHOLD 250
#define PLANE_THICKNESS_EPSILON 0.0000001f

enum PLANEINTERSECT
{
	COPLANAR,
	FRONT,
	BACK,
	STRADDLING
};

enum POINTINTERSECT
{
	POINT_IN_FRONT_OF_PLANE,
	POINT_BEHIND_PLANE,
	POINT_ON_PLANE,
};



enum BoundType
{
	Box_,
	Sphere_,
};


struct PoligonalPart
{
	Polygons* front;
	Polygons* back;
};

class BSPTree
{
public:
	BSPTree(std::vector<Polygons*> Array, int depth, int numPoly = 250);

	void Draw(Camera* cam);

private:

	void BuildBSPtree(std::vector<Polygons*>& Array, int depth, int numPoly = 250);
	Plane PickSplittingPlane(std::vector<Polygons*>& Array);

	// Return value specifying whether the polygon poly lies in front of,
	// behind of, on, or straddles the plane
	int ClassifyPolygonToPlane(Polygons* tri, Plane plane);


	// Classify point p to a plane thickened by a given thickness epsilon
	int ClassifyPointToPlane(Eigen::Vector3f p, Plane plane);

	Plane GetPlaneFromTriangle(Polygons* tri);

	void SplitTri(Polygons& tri, Plane plane, Polygons** front, Polygons** back);


	Eigen::Vector3f IntersectEdgeAgainstPlane(Eigen::Vector3f a, Eigen::Vector3f b, Plane plane);

	Node* root = nullptr;

	std::vector<Node*> leaves;

};

