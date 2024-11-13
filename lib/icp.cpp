#include "icp.h"

void scale(Mesh* mesh1, float ratio) {
    for (int v = 0; v < (int)mesh1->verts.size(); v++)
            for (int c = 0; c < 3; c++)
                mesh1->verts[v]->coords[c] *= ratio; // scaled by the euclidean ratio
}

float* calculateCenter(Mesh* mesh) {
    float *center = new float[3];
    float sum[3] = {0.0f, 0.0f, 0.0f};

    for (int i = 0; i < mesh->numVerts; i++)
        for (int c = 0; c < 3; c++)
            sum[c] += mesh->verts[i]->coords[c];

    for (int c = 0; c < 3; c++)
        center[c] = sum[c] / (float)mesh->numVerts;
    
    return center;
}

float* calculateCenter(Mesh* mesh, vector<int> indices) {
    float *center = new float[3];
    float sum[3] = {0.0f, 0.0f, 0.0f};

    for (int i = 0; i < indices.size(); i++)
        for (int c = 0; c < 3; c++)
            sum[c] += mesh->verts[indices[i]]->coords[c];

    for (int c = 0; c < 3; c++)
        center[c] = sum[c] / (float)indices.size();
    
    return center;
}

// Euclidean distance between v1 & v2
inline float distanceBetween(float *v1, float *v2)
{
	return (float)sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]) + (v2[2] - v1[2]) * (v2[2] - v1[2]));
}

// Squared Euclidean distance between v1 & v2 (no sqrt hence faster; good for comparisons only)
inline float distanceBetween2(float *v1, float *v2)
{
	return (float)(v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]) + (v2[2] - v1[2]) * (v2[2] - v1[2]);
}

// Brute force search for closest matches in quadratic time (n^2)
vector<int> searchCorrespondences(Mesh* mesh1, Mesh* mesh2, bool oneToOneMatches) {
    float dist;
    vector<int> matchedP2s;
    for (int bv1 = 0; bv1 < mesh1->numVerts; bv1++) // for each mesh verts
    {
        float minDist = INF;
        int mv = -1;
        for (int bv2 = 0; bv2 < mesh2->numVerts; bv2++) // find its closest voxel match
        {
            // omit sqrt() to save time
            dist = distanceBetween2(mesh2->verts[bv2]->coords, mesh1->verts[bv1]->coords);
            if (dist > minDist)
                continue;

            if (!oneToOneMatches || mesh2->verts[bv2]->closestVertIdx == -1)
            {
                minDist = dist;
                mv = bv2;
            }
        }
        if (mv == -1)
            cout << "WARNING: this may happen in oneToOneMatches mode if a valid mesh2 point candidate cannot be found; easy fix: oneToOneMatches=false\n";
        matchedP2s.push_back(mv);
        mesh2->verts[mv]->closestVertIdx = bv1;
    }
    return matchedP2s;
}

// transform this mesh towards the fixed mesh2 using ICP;
//   closed-form rotation matrix is from eq. 21 of the original paper: A Method for Registration of 3-D Shapes
// TRANSFORMING MESH: mesh1 -- FIXED MESH: mesh2
void ICP::align(Mesh *mesh1, Mesh *mesh2, int nMaxIters, bool oneToOne, float minDisplacement, bool prescale)
{
    cout << "Scale-Adaptive ICP in action (kd-tree not in use; affects running time drastically for large inputs)..\n";

    bool matured = false;

    // Scale to ratio between [distance between farthest 2 points on fixed mesh2] and [distance between farthest 2 points on transforming mesh1]
    if (prescale)
        scale(mesh1, (mesh2->maxEucDist / mesh1->maxEucDist)); // scaled by the euclidean ratio

    // Size of point set 1 (transforming set) and point set 2 (fixed set);
    int nP1 = mesh1->numVerts;
    int nP2 = mesh2->numVerts;

    float *m1 = new float[3];
    float *m2 = new float[3];

    MatrixXf ms1(3, nP1), ms2(3, nP1); // Mean-shifted kxnP1 matrix for the mesh1.samples and matched samples in mesh2; X is for dynamic size f is for float
    Matrix3f cov12;                    // kxk cross-covariance mtrx of 2 sets, mesh & voxels where k = 3
    Matrix4f Q;                        // complicated matrix created using cov12

    float **R1 = new float *[3]; // desired 3x3 rotation mtrx of current iteration
    for (int i = 0; i < 3; i++)
        R1[i] = new float[3];

    // Parameters of the linear system to be solved in Adaptive Scaling
    float A, B;
    float *C = new float[3];
    float *D = new float[3];

    // Final results variables (rotation, scaling and translation)
    transformations.clear();

    // Iteration parameters
    int nIters = 0;
    float dist;                                      // closest w.r.t. Euclidean distance
    float mse = INF, prevMse = 0.0f;                 // Mean-square error objective function to be minimized
    float alignErr = 100.0f, closeness = 0.0000001f; // 0.0000001f;//for face.xyz  //0.00001f;//for completehomer.xyz // 0.000000000001f default

    while (++nIters <= nMaxIters)
    {
        if (fabs(mse - prevMse) <= closeness)
        {
            // float maxDist = calculateMaxEucDist();
            // if(maxDist < mesh2->maxEucDist / 2.0) {
            // 	for (int v = 0; v < (int)verts.size(); v++)
            // 		for (int c = 0; c < 3; c++)
            // 			verts[v]->coords[c] *= (mesh2->maxEucDist / maxDist); // scaled by the euclidean ratio
            // }
            // else
            break;
        }
       
        /*
            FIND CORRESPONDENCES
        */
        // true means 1-to-1 matching will be used between this mesh points and mesh2 points,
        // i.e., no mesh2 point appears more than once in the closest-pnt matching
        // (good for scaleAdaptive since it prevents selection of the same points over and over again when the scale difference is high)
        bool oneToOneMatches = oneToOne;

        // relax one-to-one condition (= convert true to false) since the movement was so small compared to the previous iteration
        // np1 > np2 should not occur 'cos fixed mesh2 should be the bigger (more complete) one anyway
        if (fabs(mse - prevMse) < minDisplacement || matured || nP1 > nP2)
        {
            if (!matured && oneToOneMatches)
                cout << "Disabling one-to-one pairing of points (" << fabs(mse - prevMse) << ", " << nIters << ")\n";
            matured = true;
            oneToOneMatches = false;
        }

        // Get ready for the upcoming iteration; to prevent selection of a mesh2 vertex more than once
        if (oneToOneMatches)
            for (int i = 0; i < nP2; i++)
                mesh2->verts[i]->closestVertIdx = -1;

        // idxs of mesh2 points that are matched to mesh1 points
        vector<int> matchedP2s = searchCorrespondences(mesh1, mesh2, oneToOne);

        // New center of mass of the transforming mesh and fixed mesh
        m1 = calculateCenter(mesh1);
        m2 = calculateCenter(mesh2, matchedP2s);

        // Fill mean-shifted kxnP1 ms1 and ms2 matrix for the mesh1 and mesh2
        // Fill mean-shifted 3xnVerts matrix for the 2nd set (matched points, hence of size nVerts for sure)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < nP1; j++) {
                ms1(i, j) = mesh1->verts[mesh2->verts[matchedP2s[j]]->closestVertIdx]->coords[i] - m1[i];
                ms2(i, j) = mesh2->verts[matchedP2s[j]]->coords[i] - m2[i];
            }
        MatrixXf ms2t = ms2.transpose(); // transpose of ms

        // kxk cross-covariance mtrx between mesh1 and mesh2
        cov12 = ms1 * ms2t; // k x nP1 * nP1 x k
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                cov12(i, j) /= nP1;

        // Q = 4x4 complicated mtrx, see the paper
        Q(0, 0) = cov12(0, 0) + cov12(1, 1) + cov12(2, 2);
        Q(1, 0) = Q(0, 1) = cov12(1, 2) - cov12(2, 1); //(paper notation: S23 - S32, say S to cov12, and idxs of paper, i.e. start from 1)
        Q(2, 0) = Q(0, 2) = cov12(2, 0) - cov12(0, 2); // S31 - S13
        Q(3, 0) = Q(0, 3) = cov12(0, 1) - cov12(1, 0); // S12 - S21
        Q(1, 1) = 2 * cov12(0, 0) - Q(0, 0);
        Q(1, 2) = cov12(0, 1) + cov12(1, 0);
        Q(1, 3) = cov12(0, 2) + cov12(2, 0);
        Q(2, 1) = cov12(1, 0) + cov12(0, 1);
        Q(2, 2) = 2 * cov12(1, 1) - Q(0, 0);
        Q(2, 3) = cov12(1, 2) + cov12(2, 1);
        Q(3, 1) = cov12(2, 0) + cov12(0, 2);
        Q(3, 2) = cov12(2, 1) + cov12(1, 2);
        Q(3, 3) = 2 * cov12(2, 2) - Q(0, 0);

        // Desired rotation is decided by 1x4 quaternion vector, which is the (largest eigval's) eigvec of Q
        SelfAdjointEigenSolver<Matrix4f> eigensolver(Q); // solve eigenvalue/vector for symmetric n by n matrices, a.k.a. selfadjoint matrices
        if (eigensolver.info() != Success)
            abort();
        Vector4f eig_vals = eigensolver.eigenvalues();  // Given in ascending order
        Matrix4f eig_vecs = eigensolver.eigenvectors(); // Last column of eig_vecs corresponds to the largest eigenvalue, hence the one i want

        // 1x4 quaternion = [q0, q1, q2, q3]
        float q0 = eig_vecs(0, 3), q1 = eig_vecs(1, 3), q2 = eig_vecs(2, 3), q3 = eig_vecs(3, 3);
        // if (q0 < 0) { q0 = -q0;	q1 = -q1;	q2 = -q2;	q3 = -q3; } //needed??????????

        // Fill desired 3x3 rotation mtrx (i'm used to 4x4 mtrxs w/ coords but this 3x3 works because of quaternions, i guess)
        R1[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
        R1[0][1] = 2 * (q1 * q2 - q0 * q3);
        R1[0][2] = 2 * (q1 * q3 + q0 * q2);
        R1[1][0] = 2 * (q1 * q2 + q0 * q3);
        R1[1][1] = q0 * q0 + q2 * q2 - q1 * q1 - q3 * q3;
        R1[1][2] = 2 * (q2 * q3 - q0 * q1);
        R1[2][0] = 2 * (q1 * q3 - q0 * q2);
        R1[2][1] = 2 * (q2 * q3 + q0 * q1);
        R1[2][2] = q0 * q0 + q3 * q3 - q1 * q1 - q2 * q2;

        // our ICP that solves for translation and uniform scaling simulataneously to improve the initial scale guess as well

        // the math of this part is explained in my paper titled Skuller; check it out
        A = B = 0.0f;
        for (int i = 0; i < 3; i++)
            C[i] = D[i] = 0.0f;

        // Apply just ICP rotation of this iteration to all verts of the transforming mesh
        for (int v = 0; v < mesh1->numVerts; v++)
        {
            float x = mesh1->verts[v]->coords[0],
                  y = mesh1->verts[v]->coords[1],
                  z = mesh1->verts[v]->coords[2];
            // rotation
            mesh1->verts[v]->coords[0] = R1[0][0] * x + R1[0][1] * y + R1[0][2] * z;
            mesh1->verts[v]->coords[1] = R1[1][0] * x + R1[1][1] * y + R1[1][2] * z;
            mesh1->verts[v]->coords[2] = R1[2][0] * x + R1[2][1] * y + R1[2][2] * z;
        }

        for (int v = 0; v < nP1; v++)
        {
            // D = [Dx Dy Dz] is computed Dx is the sum of x-coords of all matched mesh2 points,
            //   and so on (since matched mesh2 points change each iteration D is not fixed);
            //   similarly C = [Cx Cy Cz] for transformed, i.e. ICP-rotated, mesh points
            int m2v = matchedP2s[v];
            for (int i = 0; i < 3; i++)
            {
                C[i] += mesh1->verts[v]->coords[i];
                D[i] += mesh2->verts[m2v]->coords[i];
            }
            A += mesh1->verts[v]->coords[0] * mesh1->verts[v]->coords[0] + mesh1->verts[v]->coords[1] * mesh1->verts[v]->coords[1] + mesh1->verts[v]->coords[2] * mesh1->verts[v]->coords[2];
            B += mesh1->verts[v]->coords[0] * mesh2->verts[m2v]->coords[0] + mesh1->verts[v]->coords[1] * mesh2->verts[m2v]->coords[1] + mesh1->verts[v]->coords[2] * mesh2->verts[m2v]->coords[2];
        }

        // solve for Px = r linear system to decide x = [s t1.x t1.y t1.z] where s is the uniform scale factor, t1[] is the translation
        Matrix4f P;
        P << A, C[0], C[1], C[2], C[0], (float)nP1, 0.0f, 0.0f, C[1], 0.0f, (float)nP1, 0.0f, C[2], 0.0f, 0.0f, (float)nP1;
        // cout << "Here is the matrix P:\n" << P << endl;
        Vector4f r;
        r << B, D[0], D[1], D[2];
        // cout << "Here is the right hand side r:\n" << r << endl;
        Vector4f x = P.jacobiSvd(ComputeFullU | ComputeFullV).solve(r); // s = x(0) is the uniform scale factor, and the rest is the translation t1
        // cout << "The least-squares solution is:\n" << x << "\nsame as r?\n" << P * x << endl;

        // apply translation and scaling to all vertices
        for (int v = 0; v < mesh1->numVerts; v++)
        {
            // scaling
            mesh1->verts[v]->coords[0] *= x(0);
            mesh1->verts[v]->coords[1] *= x(0);
            mesh1->verts[v]->coords[2] *= x(0);

            // translation
            mesh1->verts[v]->coords[0] += x(1);
            mesh1->verts[v]->coords[1] += x(2);
            mesh1->verts[v]->coords[2] += x(3);
        }
        // cout << x(0) << "\t\t" << x(1) << " " << x(2) << " " << x(3) << "\n";//cout << x(0) << " ";

        // Update final values of transformation
        Vector4f xIt = x;
        float **RIt = new float *[3]; // desired 3x3 rotation mtrx of current iteration
        for (int i = 0; i < 3; i++)
        {
            RIt[i] = new float[3];
            for (int j = 0; j < 3; j++)
                RIt[i][j] = R1[i][j];
        }
        transformations.push_back(make_pair(xIt, RIt));

        prevMse = mse; // mse of previous iteration as a termination condition
        // compute new mse value as an iteration condition
        mse = 0.0f;
        for (int bv2 = 0; bv2 < nP1; bv2++)
        {
            // squared distance b/w (target) matched voxel and its closest match in mesh whose coords just updated above
            dist = pow(mesh2->verts[matchedP2s[bv2]]->coords[0] - mesh1->verts[mesh2->verts[matchedP2s[bv2]]->closestVertIdx]->coords[0], 2.0f) +
                   pow(mesh2->verts[matchedP2s[bv2]]->coords[1] - mesh1->verts[mesh2->verts[matchedP2s[bv2]]->closestVertIdx]->coords[1], 2.0f) +
                   pow(mesh2->verts[matchedP2s[bv2]]->coords[2] - mesh1->verts[mesh2->verts[matchedP2s[bv2]]->closestVertIdx]->coords[2], 2.0f);
            mse += dist; // no sqrt() above for efficiency, hence dist = squared dist indeed
        }
        mse /= nP1;
        if (nIters % 10 == 1)
            cout << "mse of ICP iter" << nIters - 1 << ": " << mse << "\tmse diff: " << fabs(mse - prevMse) << endl;

    } // end of while nIters

    // recapture memo
    for (int i = 0; i < 3; i++)
        delete[] R1[i];
    delete[] R1;
    delete[] C;
    delete[] D;

    cout << nIters - 1 << "'th ICP iteration w/ final mse: " << mse << "\nICP done!\n";
}

void ICP::transform(Mesh* mesh)
{
    for(int i=0; i < transformations.size(); i++) {
        Vector4f st = transformations[i].first;
        float ** r = transformations[i].second;

        for (int v = 0; v < mesh->numVerts; v++)
        {
            float x = mesh->verts[v]->coords[0],
                y = mesh->verts[v]->coords[1],
                z = mesh->verts[v]->coords[2];
            // rotation
            mesh->verts[v]->coords[0] = r[0][0] * x + r[0][1] * y + r[0][2] * z;
            mesh->verts[v]->coords[1] = r[1][0] * x + r[1][1] * y + r[1][2] * z;
            mesh->verts[v]->coords[2] = r[2][0] * x + r[2][1] * y + r[2][2] * z;
        }
        for (int v = 0; v < mesh->numVerts; v++)
        {
            float x = mesh->verts[v]->coords[0],
                y = mesh->verts[v]->coords[1],
                z = mesh->verts[v]->coords[2];
            // scaling
            mesh->verts[v]->coords[0] *= st(0);
            mesh->verts[v]->coords[1] *= st(0);
            mesh->verts[v]->coords[2] *= st(0);
            // translation
            mesh->verts[v]->coords[0] += st(1);
            mesh->verts[v]->coords[1] += st(2);
            mesh->verts[v]->coords[2] += st(3);
        }
    }
}