#pragma once
#include "BoundingVolume.h"

class RitterSphere : public BoundingVolume
{
public:
	RitterSphere(std::vector<Eigen::Vector3f>& verts, std::string file = "../common/models/boundingvolumes/sphere.obj");
private:
	void MostSeparatedPointsOnAABB(int& min, int& max, std::vector<Eigen::Vector3f>& verts, std::vector<Eigen::Vector3f>& ogVerts);
	void SphereFromDistantPoints(Eigen::Vector3f& center, float& rad, std::vector<Eigen::Vector3f>& verts, std::vector<Eigen::Vector3f>& ogVerts);
	void CalcRitterSphere(Eigen::Vector3f& center, float& rad, std::vector<Eigen::Vector3f>& verts, std::vector<Eigen::Vector3f>& ogverts);
	void RitterIterative(Eigen::Vector3f& center, float& rad, std::vector<Eigen::Vector3f>& ogVerts, std::vector<Eigen::Vector3f>& verts_);
};