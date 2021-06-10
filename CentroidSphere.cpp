#include "CentroidSphere.h"

CentroidSphere::CentroidSphere(std::vector<Eigen::Vector3f>& verts, std::string file): BoundingVolume(file)
{
    float INF = std::numeric_limits<float>::infinity();
    Eigen::Vector3f minV(INF, INF, INF);
    Eigen::Vector3f maxV(-INF, -INF, -INF);
    Eigen::Vector3f centroid(0, 0, 0);

    for (unsigned i = 0; i < verts.size(); ++i)
    {
        centroid += verts[i];
        minV.x() = std::min(minV.x(), verts[i].x());
        minV.y() = std::min(minV.y(), verts[i].y());
        minV.z() = std::min(minV.z(), verts[i].z());

        maxV.x() = std::max(maxV.x(), verts[i].x());
        maxV.y() = std::max(maxV.y(), verts[i].y());
        maxV.z() = std::max(maxV.z(), verts[i].z());
    }
    centroid /= verts.size();
    float r = std::max((centroid - minV).norm(), (centroid - maxV).norm());

    Transform* trans = reinterpret_cast<Transform*>(components[Transform_]);


    Eigen::Vector3f scale = Eigen::Vector3f(0, 0, r * 1.2);
    trans->Scale(r * 1.2);
    trans->setPosition(((maxV + minV) / 2.f) + scale);
    trans->Update();
}

CentroidSphere::CentroidSphere(Eigen::Vector3f min, Eigen::Vector3f max, std::string file) : BoundingVolume(file)
{
    float INF = std::numeric_limits<float>::infinity();
    Eigen::Vector3f minV = min;
    Eigen::Vector3f maxV = max;
    Eigen::Vector3f centroid = (maxV + minV ) / 2.f;
    float r = std::max((centroid - minV).norm(), (centroid - maxV).norm());

    Transform* trans = reinterpret_cast<Transform*>(components[Transform_]);

    Eigen::Vector3f scale = Eigen::Vector3f(0, 0, r) ;
    trans->Scale(r);
    trans->setPosition(centroid + scale);
    trans->Update();
}
