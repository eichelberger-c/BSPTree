#ifndef BOUNDING_VOLUME_HEADER
#define BOUNDING_VOLUME_HEADER

#include "stdafx.h"
#include <GL/glew.h>
#include "ShaderManager.h"
#include "Camera.h"
#include "Transform.h"
#include "Object.h"
enum BoundingTypes
{
    AABB_,
    Centroid_,
    PCA_,
    Ritter_,
    OBB_,
    BVNumTypes_
};

class BoundingVolume : public Obj
{
public:
    BoundingVolume(std::string name) : Obj(name, false) {}
    Eigen::Vector3f& GetCenter() { return m_center; }
protected:

    Eigen::Vector3f m_center;
};

// Given Sphere s and Point p, update s (if needed) to just encompass p
void SphereOfSphereAndPt(Eigen::Vector3f& center, float& rad, Eigen::Vector3f& p);
#endif