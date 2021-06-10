#include "BSPTree.h"

BSPTree::BSPTree(std::vector<Polygons*> Array, int depth, int numPoly)
{
    BuildBSPtree(Array, depth, numPoly);
}

void BSPTree::Draw(Camera* cam)
{
    for (int i = 0; i < leaves.size(); ++i)
    {
        leaves[i]->Draw(cam);
    }
}

void BSPTree::BuildBSPtree(std::vector<Polygons*>& polygons, int depth, int numPoly)
{
    // Return NULL tree if there are no polygons
    if (polygons.empty()) return;

    // Get number of polygons in the input vector
    int numPolygons = polygons.size();

    // If criterion for a leaf is matched, create a leaf node from remaining polygons
    if (numPolygons <= numPoly)
    {
        Node* leaf = new Node(polygons);
        leaves.push_back(leaf);
        return;
        //return leaf;
    }

    // Select best possible partitioning plane based on the input geometry
    Plane splitPlane = PickSplittingPlane(polygons);

    std::vector<Polygons*> frontList, backList;

    // Test each polygon against the dividing plane, adding them
    // to the front list, back list, or both, as appropriate
    for (int i = 0; i < numPolygons; i++)
    {
        Polygons* poly = polygons[i];
        Polygons* front, *back;
        switch (ClassifyPolygonToPlane(poly, splitPlane))
        {
        case COPLANAR:
        case FRONT:
            frontList.push_back(poly);
            break;
        case BACK:
            backList.push_back(poly);
            break;
        case STRADDLING:
            //// Split polygon to plane and send a part to each side of the plane
            SplitTri(*poly, splitPlane, &front, &back);
            frontList.push_back(front);
            backList.push_back(back);
            break;
        }
    }
    BuildBSPtree(frontList, depth + 1, numPoly); 
    BuildBSPtree(backList, depth + 1, numPoly);
    //return new Node(frontTree, backTree);
}

Plane BSPTree::PickSplittingPlane(std::vector<Polygons*>& Array)
{
    // Blend factor for optimizing for balance or splits (should be tweaked)
    const float K = .8f;
    // Variables for tracking best splitting plane seen so far
    Plane bestPlane;
    float bestScore = FLT_MAX;
    // Try the plane of each polygon as a dividing plane
    for (int i = 0; i < Array.size(); i ++)
    {
        int numInFront = 0, numBehind = 0, numStraddling = 0;
        Plane plane = GetPlaneFromTriangle(Array[i]);
        // Test against all other polygons
        for (int j = 0; j < Array.size(); j++) {
            // Ignore testing against self
            if (i == j) continue;
            // Keep standing count of the various poly-plane relationships
            switch (ClassifyPolygonToPlane(Array[j], plane))
            {
            case COPLANAR:
                /* Coplanar polygons treated as being in front of plane */
            case FRONT:
                numInFront++;
                break;
            case BACK:
                numBehind++;
                break;
            case STRADDLING:
                numStraddling++;
                break;
            }
        }
        // Compute score as a weighted combination (based on K, with K in range
        // 0..1) between balance and splits (lower score is better)
        float score = K * numStraddling + (1.0f - K) * abs(numInFront - numBehind);
        if (score < bestScore)
        {
            bestScore = score;
            bestPlane = plane;
        }
    }
    return bestPlane;
}

// Return value specifying whether the polygon poly lies in front of,
// behind of, on, or straddles the plane
int BSPTree::ClassifyPolygonToPlane(Polygons* tri, Plane plane)
{
    // Loop over all polygon vertices and count how many vertices
    // lie in front of and how many lie behind of the thickened plane
    int numInFront = 0;
    int numBehind = 0;
    int numVerts = 1;
    for (int i = 0; i < numVerts; i++)
    {
        Eigen::Vector3f p = tri->GetVertex(i);
        switch (ClassifyPointToPlane(p, plane)) {
        case POINT_IN_FRONT_OF_PLANE:
            numInFront++;
            break;
        case POINT_BEHIND_PLANE:
            numBehind++;
            break;
        }
    }

    // If vertices on both sides of the plane, the polygon is straddling
    if (numBehind != 0 && numInFront != 0)
        return STRADDLING;
    // If one or more vertices in front of the plane and no vertices behind
    // the plane, the polygon lies in front of the plane
    if (numInFront != 0)
        return FRONT;
    // Ditto, the polygon lies behind the plane if no vertices in front of
    // the plane, and one or more vertices behind the plane
    if (numBehind != 0)
        return BACK;
    // All vertices lie on the plane so the polygon is coplanar with the plane
    return COPLANAR;
}


// Classify point p to a plane thickened by a given thickness epsilon
int BSPTree::ClassifyPointToPlane(Eigen::Vector3f p, Plane plane)
{
    float dist = plane.n.dot(p) - plane.d;
    // Classify p based on the signed distance
    if (dist > 0.0000001f)
        return POINT_IN_FRONT_OF_PLANE;
    if (dist < 0.0000001f)
        return POINT_BEHIND_PLANE;
    return POINT_ON_PLANE;
}


Plane BSPTree::GetPlaneFromTriangle(Polygons* tri)
{
    Eigen::Vector3f U;
    U = tri->GetVertex(1) - tri->GetVertex(0);

    Eigen::Vector3f V;
    V = tri->GetVertex(2) - tri->GetVertex(0);

    Eigen::Vector3f N;
    N = U.cross(V).normalized();

    Plane p;
    p.n = N;

    p.pos = tri->GetVertex(0);

    p.d = (p.pos.x() * N.x() + p.pos.y() * N.y() + p.pos.z() * N.z());

    return p;
}


void BSPTree::SplitTri(Polygons& poly, Plane plane, Polygons** front, Polygons** back)
{
    int numFront = 0, numBack = 0;
    std::vector<Eigen::Vector3f> frontVerts, backVerts;
    // Test all edges (a, b) starting with edge from last to first vertex
    int numVerts = poly.m_verts.size();
    Eigen::Vector3f a = poly.GetVertex(poly.m_verts.size()-1);
    int aSide = ClassifyPointToPlane(a, plane);

    // Loop over all edges given by vertex pair (n - 1, n)
    for (int n = 0; n < numVerts; n++)
    {
        Eigen::Vector3f b = poly.GetVertex(n);
        int bSide = ClassifyPointToPlane(b, plane);
        if (bSide == POINT_IN_FRONT_OF_PLANE)
        {
            if (aSide == POINT_BEHIND_PLANE)
            {
                // Edge (a, b) straddles, output intersection point to both sides
                Eigen::Vector3f i = IntersectEdgeAgainstPlane(b, a, plane);
                (frontVerts).push_back(i);
                (backVerts).push_back(i);
            }
            // In all three cases, output b to the front side
            (frontVerts).push_back(b);
        }
        else if (bSide == POINT_BEHIND_PLANE)
        {
            if (aSide == POINT_IN_FRONT_OF_PLANE)
            {
                // Edge (a, b) straddles plane, output intersection point
                Eigen::Vector3f i = IntersectEdgeAgainstPlane(a, b, plane);
                (frontVerts).push_back(i);
                (backVerts).push_back(i);
            }
            else if (aSide == POINT_ON_PLANE)
            {
                // Output a when edge (a, b) goes from ‘on’ to ‘behind’ plane
                (backVerts).push_back(a);
            }
            // In all three cases, output b to the back side
            (backVerts).push_back(b);
        }
        else
        {
            // b is on the plane. In all three cases output b to the front side
            (frontVerts).push_back(b);
            // In one case, also output b to back side
            if (aSide == POINT_BEHIND_PLANE)
                (backVerts).push_back(b);
        }
        // Keep b as the starting point of the next edge
        a = b;
        aSide = bSide;
    }
    *front = new Polygons(frontVerts); 
    *back = new Polygons(backVerts);
}


Eigen::Vector3f BSPTree::IntersectEdgeAgainstPlane(Eigen::Vector3f a, Eigen::Vector3f b, Plane plane)
{
    Eigen::Vector3f dir = b - a;

    float first = (plane.pos - a).dot(plane.n);
    float third = dir.dot(plane.n);

    float t = (first) / third;

    return a + dir * t;
}