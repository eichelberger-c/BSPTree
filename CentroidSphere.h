#pragma once
#include "BoundingVolume.h"

class CentroidSphere : public BoundingVolume
{
public:
	CentroidSphere(std::vector<Eigen::Vector3f>& verts, std::string file = "../common/models/boundingvolumes/sphere.obj");
	CentroidSphere(Eigen::Vector3f min, Eigen::Vector3f max, std::string file = "../common/models/boundingvolumes/sphere.obj");
private:
};