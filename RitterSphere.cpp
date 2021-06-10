#include "RitterSphere.h"

RitterSphere::RitterSphere(std::vector<Eigen::Vector3f>& verts, std::string file) : BoundingVolume(file)
{
    Eigen::Vector3f center;
    float rad;
    std::vector<Eigen::Vector3f> tempVerts = verts;
	RitterIterative(center, rad, verts, tempVerts);

    Transform* trans = reinterpret_cast<Transform*>(components[Transform_]);

    Eigen::Vector3f scale = Eigen::Vector3f(0, 0, rad * 1.2);
    trans->Scale(rad * 1.2);
    trans->setPosition(center + scale);
    trans->Update();
}

int random(int min, int max) //range : [min, max)
{
    static bool first = true;
    if (first)
    {
        srand(time(NULL)); //seeding for the first time only!
        first = false;
    }
    int ret = min + rand() % ((max + 1) - min);

    if (ret >= max)
        return max - 1;
    if (ret <= min)
        return min;

    return ret;
}


// Compute indices to the two most separated points of the (up to) six points
// defining the AABB encompassing the point set. Return these as min and max.
void RitterSphere::MostSeparatedPointsOnAABB(int& min, int& max, std::vector<Eigen::Vector3f>& verts, std::vector<Eigen::Vector3f>& ogVerts)
{
    // First find most extreme points along principal axes
    int minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
    for (int i = 1; i < ogVerts.size(); i++) {
        if (verts[i].x() < verts[minx].x()) minx = i;
        if (verts[i].x() > verts[maxx].x()) maxx = i;
        if (verts[i].y() < verts[miny].y()) miny = i;
        if (verts[i].y() > verts[maxy].y()) maxy = i;
        if (verts[i].z() < verts[minz].z()) minz = i;
        if (verts[i].z() > verts[maxz].z()) maxz = i;
    }
    // Compute the squared distances for the three pairs of points
    float dist2x = (verts[maxx] - verts[minx]).dot(verts[maxx] - verts[minx]);
    float dist2y = (verts[maxy] - verts[miny]).dot(verts[maxy] - verts[miny]);
    float dist2z = (verts[maxz] - verts[minz]).dot(verts[maxz] - verts[minz]);
    // Pick the pair (min,max) of points most distant
    min = minx;
    max = maxx;
    if (dist2y > dist2x&& dist2y > dist2z) {
        max = maxy;
        min = miny;
    }
    if (dist2z > dist2x&& dist2z > dist2y) {
        max = maxz;
        min = minz;
    }
}

void RitterSphere::SphereFromDistantPoints(Eigen::Vector3f& center, float& rad, std::vector<Eigen::Vector3f>& verts, std::vector<Eigen::Vector3f>& ogVerts)
{
    // Find the most separated point pair defining the encompassing AABB
    int min, max;
    MostSeparatedPointsOnAABB(min, max, verts, ogVerts);
    // Set up sphere to just encompass these two points
    center = (verts[min] + verts[max]) * 0.5f;
    rad = (verts[max] - center).dot(verts[max] - center);
    rad = sqrt(rad);
}

void RitterSphere::CalcRitterSphere(Eigen::Vector3f& center, float& rad, std::vector<Eigen::Vector3f>& verts, std::vector<Eigen::Vector3f>& ogverts)
{
    // Get sphere encompassing two approximately most distant points
    SphereFromDistantPoints(center, rad, verts, ogverts);
    // Grow sphere to include all points
    for (int i = 0; i < verts.size(); i++)
        SphereOfSphereAndPt(center, rad, verts[i]);
}

void RitterSphere::RitterIterative(Eigen::Vector3f& center, float& rad, std::vector<Eigen::Vector3f>& ogVerts, std::vector<Eigen::Vector3f>& verts_)
{
    const int NUM_ITER = 8;
    CalcRitterSphere(center, rad, ogVerts, ogVerts);
    float rad2 = 0.f;
    Eigen::Vector3f center2;
    for (int k = 0; k < NUM_ITER; k++) {
        // Shrink sphere somewhat to make it an underestimate (not bound)
        rad2 = rad2 * 0.95f;
        // Make sphere bound data again
        for (int i = 0; i < verts_.size(); i++) {
            // Swap pt[i] with pt[j], where j randomly from interval [i+1,numPts-1]
            int j = random(i + 1, ogVerts.size());
            Eigen::Vector3f temp = verts_[i];
            verts_[i] = verts_[j];
            verts_[j] = temp;
            SphereOfSphereAndPt(center2, rad2, ogVerts[i]);
        }
        // Update s whenever a tighter sphere is found
        if (rad2 < rad)
        {
            rad = rad2;
            center = center2;
        }
    }
}
