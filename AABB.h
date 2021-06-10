#pragma once
#include "BoundingVolume.h"
class AABB : public BoundingVolume 
{
public:
	AABB(std::vector<Eigen::Vector3f>& verts, std::string name = "../common/models/boundingvolumes/cube.obj");
	AABB(Eigen::Vector3f min, Eigen::Vector3f max, std::string name = "../common/models/boundingvolumes/cube.obj");
};