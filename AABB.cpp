#include "AABB.h"

#include "Model.h"
AABB::AABB(std::vector<Eigen::Vector3f>& verts, std::string name) : BoundingVolume(name)
{
    Model* model = reinterpret_cast<Model*>(components[Model_]);
    //model->SetMode(GL_LINES);
    Transform* trans = reinterpret_cast<Transform*>(components[Transform_]);
    //m_modelVerts = verts;
    float INF = std::numeric_limits<float>::infinity();
    Eigen::Vector3f minV(INF, INF, INF);
    Eigen::Vector3f maxV(-INF, -INF, -INF);
    for (unsigned i = 0; i < verts.size(); ++i)
    {
        minV.x() = std::min(minV.x(), verts[i].x());
        minV.y() = std::min(minV.y(), verts[i].y());
        minV.z() = std::min(minV.z(), verts[i].z());

        maxV.x() = std::max(maxV.x(), verts[i].x());
        maxV.y() = std::max(maxV.y(), verts[i].y());
        maxV.z() = std::max(maxV.z(), verts[i].z());

    }
    float width, height, depth;

    height = std::abs(maxV.y() - minV.y());
    width = std::abs(maxV.x() - minV.x());
    depth = std::abs(maxV.z() - minV.z());
    Eigen::Vector3f center = (maxV + minV) / 2.f;
    //Eigen::Vector3f center = (maxV + minV) / 2.f;
    Eigen::Vector3f scale = Eigen::Vector3f(width, height, depth) / 2.f;
    trans->setPosition(center - scale);
    trans->setSclae(scale);
    trans->Update();
}

AABB::AABB(Eigen::Vector3f min, Eigen::Vector3f max, std::string name) : BoundingVolume(name)
{
    Model* model = reinterpret_cast<Model*>(components[Model_]);
    //model->SetMode(GL_LINES);
    Transform* trans = reinterpret_cast<Transform*>(components[Transform_]);
    //m_modelVerts = verts;
    float width, height, depth;

    height = std::abs(max.y() - min.y());
    width = std::abs(max.x() - min.x());
    depth = std::abs(max.z() - min.z());

    Eigen::Vector3f center = (max + min) / 2.f;
    Eigen::Vector3f scale = Eigen::Vector3f(width, height, depth) / 2.f;
    trans->setPosition(center - scale);
    trans->setSclae(scale);
    trans->Update();
}