#include "BoundingVolume.h"

void CalcSphere(float rad, Eigen::Vector3f center, std::vector<Eigen::Vector3f>& vertices, std::vector<unsigned>& faces)
{
    float x, y, z, xy;                              // vertex position
    float numDivisions2 = float(360) / 32;
    float sectorCount = float(180) / 32;
    float stackCount = float(360) / 32;

    float sectorStep = 2 * M_PI / sectorCount;
    float stackStep = M_PI / stackCount;
    float sectorAngle, stackAngle;

    for (int i = 0; i <= stackCount; ++i)
    {
        stackAngle = M_PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = rad * cosf(stackAngle);             // r * cos(u)
        z = rad * sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for (int j = 0; j <= sectorCount; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position (x, y, z)
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            vertices.push_back(Eigen::Vector3f(x, y, z) + center);
        }
    }


    int k1, k2;
    for (int i = 0; i < stackCount; ++i)
    {
        k1 = int(i * (sectorCount + 1));     // beginning of current stack
        k2 = int(k1 + sectorCount + 1);      // beginning of next stack

        for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding first and last stacks
            // k1 => k2 => k1+1
            if (i != 0)
            {
                faces.push_back(k1);
                faces.push_back(k2);
                faces.push_back(k1 + 1);
            }

            // k1+1 => k2 => k2+1
            if (i != (stackCount - 1))
            {
                faces.push_back(k1 + 1);
                faces.push_back(k2);
                faces.push_back(k2 + 1);
            }
        }
    }
}

void SphereOfSphereAndPt(Eigen::Vector3f& center, float& rad, Eigen::Vector3f& p)
{
    // Compute squared distance between point and sphere center
    Eigen::Vector3f d = p - center;
    float dist2 = d.dot(d);
    // Only update s if point p is outside it
    if (dist2 > rad* rad) {
        float dist = sqrt(dist2);
        float newRadius = (rad + dist) * 0.5f;
        float k = (newRadius - rad) / dist;
        rad = newRadius;
        center += d * k;
    }
}
