#include "PCASphere.h"

PCASphere::PCASphere(std::vector<Eigen::Vector3f>& verts, std::string file) : BoundingVolume(file)
{
    Eigen::Vector3f center;
    float rad;
	RitterEigenSphere(verts, center, rad);

    Transform* trans = reinterpret_cast<Transform*>(components[Transform_]);

    Eigen::Vector3f scale = Eigen::Vector3f(0, 0, rad * 1.2);
    trans->Scale(rad * 1.2);
    trans->setPosition(center + scale);
    trans->Update();
}


// Compute variance of a set of 1D values
float PCASphere::Variance(float x[], int n)
{
    float u = 0.0f;
    for (int i = 0; i < n; i++)
        u += x[i];
    u /= n;
    float s2 = 0.0f;
    for (int i = 0; i < n; i++)
        s2 += (x[i] - u) * (x[i] - u);
    return s2 / n;
}

void PCASphere::CovarianceMatrix(Eigen::Matrix3f& cov, std::vector<Eigen::Vector3f>& verts)
{
    float oon = 1.0f / (float)verts.size();
    Eigen::Vector3f c = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
    float e00, e11, e22, e01, e02, e12;
    // Compute the center of mass (centroid) of the points
    for (int i = 0; i < verts.size(); i++)
        c += verts[i];
    c *= oon;
    // Compute covariance elements
    e00 = e11 = e22 = e01 = e02 = e12 = 0.0f;
    for (int i = 0; i < verts.size(); i++) {
        // Translate points so center of mass is at origin
        Eigen::Vector3f p = verts[i] - c;
        // Compute covariance of translated points
        e00 += p.x() * p.x();
        e11 += p.y() * p.y();
        e22 += p.z() * p.z();
        e01 += p.x() * p.y();
        e02 += p.x() * p.z();
        e12 += p.y() * p.z();
    }
    // Fill in the covariance matrix elements
    cov(0, 0) = e00 * oon;
    cov(1, 1) = e11 * oon;
    cov(2, 2) = e22 * oon;
    cov(0, 1) = cov(1, 0) = e01 * oon;
    cov(0, 2) = cov(2, 0) = e02 * oon;
    cov(1, 2) = cov(2, 1) = e12 * oon;
}


// 2-by-2 Symmetric Schur decomposition. Given an n-by-n symmetric matrix
// and indices p, q such that 1 <= p < q <= n, computes a sine-cosine pair
// (s, c) that will serve to form a Jacobi rotation matrix.
//
// See Golub, Van Loan, Matrix Computations, 3rd ed, p428
void PCASphere::SymSchur2(Eigen::Matrix3f& a, int p, int q, float& c, float& s)
{
    if (abs(a(p, q)) > 0.0001f) {
        float r = (a(q, q) - a(p, p)) / (2.0f * a(p, q));
        float t;
        if (r >= 0.0f)
            t = 1.0f / (r + sqrt(1.0f + r * r));
        else
            t = -1.0f / (-r + sqrt(1.0f + r * r));
        c = 1.0f / sqrt(1.0f + t * t);
        s = t * c;
    }
    else {
        c = 1.0f;
        s = 0.0f;
    }
}

// Computes the eigenvectors and eigenvalues of the symmetric matrix A using
// the classic Jacobi method of iteratively updating A asA=J∧T * A * J,
// where J = J(p, q, theta) is the Jacobi rotation matrix.
//
// On exit, v will contain the eigenvectors, and the diagonal elements
// of a are the corresponding eigenvalues.
//
// See Golub, Van Loan, Matrix Computations, 3rd ed, p428
void PCASphere::Jacobi(Eigen::Matrix3f& a, Eigen::Matrix3f& v)
{
    int i, j, n, p, q;
    float prevoff, c, s;
    Eigen::Matrix3f J, b, t;
    // Initialize v to identify matrix
    for (i = 0; i < 3; i++) {
        v(i, 0) = v(i, 1) = v(i, 2) = 0.0f;
        v(i, i) = 1.0f;
    }
    // Repeat for some maximum number of iterations
    const int MAX_ITERATIONS = 50;
    for (n = 0; n < MAX_ITERATIONS; n++) {
        // Find largest off-diagonal absolute element a[p][q]
        p = 0; q = 1;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i == j) continue;
                if (abs(a(i, j)) > abs(a(p, q))) {
                    p = i;
                    q = j;
                }
            }
        }
        // Compute the Jacobi rotation matrix J(p, q, theta)
        // (This code can be optimized for the three different cases of rotation)
        SymSchur2(a, p, q, c, s);
        for (i = 0; i < 3; i++) {
            J(i, 0) = J(i, 1) = J(i, 2) = 0.0f;
            J(i, i) = 1.0f;
        }
        J(p, p) = c; J(p, q) = s;
        J(q, p) = -s; J(q, q) = c;
        // Cumulate rotations into what will contain the eigenvectors
        v = v * J;
        // Make ’a’ more diagonal, until just eigenvalues remain on diagonal
        a = (J.transpose() * a) * J;
        // Compute "norm" of off-diagonal elements
        float off = 0.0f;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i == j) continue;
                off += a(i, j) * a(i, j);
            }
        }
        /* off = sqrt(off); not needed for norm comparison */
        // Stop when norm no longer decreasing
        if (n > 2 && off >= prevoff)
            return;
        prevoff = off;
    }
}

// Returns indices imin and imax into pt[] array of the least and
// most, respectively, distant points along the direction dir
void PCASphere::ExtremePointsAlongDirection(Eigen::Vector3f dir, int* imin, int* imax, std::vector<Eigen::Vector3f>& verts)
{
    float minproj = FLT_MAX, maxproj = -FLT_MAX;
    for (int i = 0; i < verts.size(); i++) {
        // Project vector from origin to point onto direction vector
        float proj = verts[i].dot(dir);
        // Keep track of least distant point along direction vector
        if (proj < minproj) {
            minproj = proj;
            *imin = i;
        }
        // Keep track of most distant point along direction vector
        if (proj > maxproj) {
            maxproj = proj;
            *imax = i;
        }
    }
}

void PCASphere::EigenSphere(std::vector<Eigen::Vector3f>& verts, Eigen::Vector3f& center, float& rad)
{
    Eigen::Matrix3f m, v;
    // Compute the covariance matrix m
    CovarianceMatrix(m, verts);
    // Decompose it into eigenvectors (in v) and eigenvalues (in m)
    Jacobi(m, v);
    // Find the component with largest magnitude eigenvalue (largest spread)
    Eigen::Vector3f e;
    int maxc = 0;
    float maxf, maxe = abs(m(0, 0));
    if ((maxf = abs(m(1, 1))) > maxe) maxc = 1, maxe = maxf;
    if ((maxf = abs(m(2, 2))) > maxe) maxc = 2, maxe = maxf;
    e[0] = v(0, maxc);
    e[1] = v(1, maxc);
    e[2] = v(2, maxc);

    int imin, imax;
    ExtremePointsAlongDirection(e, &imin, &imax, verts);
    Eigen::Vector3f minpt = verts[imin];
    Eigen::Vector3f maxpt = verts[imax];
    float dist = sqrt((maxpt - minpt).dot(maxpt - minpt));
    rad = dist * 0.5f;
    center = (minpt + maxpt) * 0.5f;
}

void PCASphere::RitterEigenSphere(std::vector<Eigen::Vector3f>& verts, Eigen::Vector3f& center, float& rad)
{
    // Start with sphere from maximum spread
    EigenSphere(verts, center, rad);
    // Grow sphere to include all points
    for (int i = 0; i < verts.size(); i++)
        SphereOfSphereAndPt(center, rad, verts[i]);
}
