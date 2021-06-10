#pragma once
#include "BoundingVolume.h"

class PCASphere : public BoundingVolume
{
public:
	PCASphere(std::vector<Eigen::Vector3f>& verts, std::string file = "../common/models/boundingvolumes/sphere.obj");
private:
	float Variance(float x[], int n);
	void CovarianceMatrix(Eigen::Matrix3f& cov, std::vector<Eigen::Vector3f>& verts);
	void SymSchur2(Eigen::Matrix3f& a, int p, int q, float& c, float& s);
	void Jacobi(Eigen::Matrix3f& a, Eigen::Matrix3f& v);
	void ExtremePointsAlongDirection(Eigen::Vector3f dir, int* imin, int* imax, std::vector<Eigen::Vector3f>& verts);
	void EigenSphere(std::vector<Eigen::Vector3f>& verts, Eigen::Vector3f& center, float& rad);
	void RitterEigenSphere(std::vector<Eigen::Vector3f>& verts, Eigen::Vector3f& center, float& rad);
};