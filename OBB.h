#pragma once
#include "BoundingVolume.h"
#include "stdafx.h"
class OBB : public BoundingVolume
{
public:
	OBB(std::vector<Eigen::Vector3f>& verts, std::string name = "../common/models/boundingvolumes/cube.obj");
private:
	void CovarianceMatrix(Eigen::Matrix3f& cov, std::vector<Eigen::Vector3f>& verts);
	void SymSchur2(Eigen::Matrix3f& a, int p, int q, float& c, float& s);
	void Jacobi(Eigen::Matrix3f& a, Eigen::Matrix3f& v);
};