#include "Mesh.h"
#include <cmath>

//Euclidean distance between v1 & v2
inline float distanceBetween(float* v1, float* v2)
{
	return (float) sqrt( (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}

// Squared Euclidean distance between v1 & v2 (no sqrt hence faster; good for comparisons only)
inline float distanceBetween2(float* v1, float* v2)
{
	return (float) (v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]);
}

bool Mesh::loadPnt(char *meshFile)
{
    cout << "Point Cloud initializing (to " << meshFile << ")... \n";

    FILE* fPtr;
    if (! (fPtr = fopen(meshFile, "r")))
	{
		cout << "cannot read " << meshFile << endl
			 << "WARNING: all input files must be in a folder named \"input\"; similarly there must exist a folder named \"output\" to write the results to\n";
		exit(0);
	}

    while (! feof( fPtr ) )
	{
		float a, b, c;
		fscanf(fPtr, "%f %f %f\n", &a, &b, &c);
		float* co = new float[3];
		co[0] = a;
		co[1] = b;
		co[2] = c;

		addVertex(co); // No duplicate check
	}

    // Can be improved (max(maxEucDist, distanceBetween(v1, v2)))
    maxEucDist = 0.0f;
	for (int i = 0; i < (int) verts.size(); i++)
		for (int j = 0; j < (int) verts.size(); j++)
			if (distanceBetween(verts[i]->coords, verts[j]->coords) > maxEucDist)
				maxEucDist = distanceBetween(verts[i]->coords, verts[j]->coords);
	cout << "distance between farthest 2 points: " << maxEucDist << endl;

	cout << "Point cloud has " << (int) verts.size() << " verts\n\n";
	return true;
}

// Fastly add verts w/ coords c to mesh; No duplication check
void Mesh::addVertex(float* c)
{
	int vSize = (int) verts.size(); //size before push_back just-below
	verts.push_back(new Vertex(vSize, c));
}

// transform this mesh towards the fixed mesh2 using ICP;
//   closed-form rotation matrix is from eq. 21 of the original paper: A Method for Registration of 3-D Shapes
// TRANSFORMING MESH: mesh1 -- FIXED MESH: mesh2
float Mesh::ICP(Mesh* mesh2, int nMaxIters, bool oneToOne, float minDisplacement)
{
    cout << "Scale-Adaptive ICP in action (kd-tree not in use; affects running time drastically for large inputs)..\n";

    bool matured = false;
    bool prescale = false; // true to prescale input pairs

	// Scale to ratio between [distance between farthest 2 points on fixed mesh2] and [distance between farthest 2 points on transforming mesh1]
    if (prescale)		
		for (int v = 0; v < (int) verts.size(); v++)
			for (int c = 0; c < 3; c++)
				verts[v]->coords[c] *= (mesh2->maxEucDist / maxEucDist); // scaled by the euclidean ratio

    // Put all verts in samples[] because code below uses samples[] everywhere.
    //  protect original samples[] in copy vectors
    for (int i = 0; i < (int) verts.size(); i++)
        samples.push_back(i);
    mesh2->samples.clear();
    for (int i = 0; i < (int) mesh2->verts.size(); i++)
        mesh2->samples.push_back(i);

    // Center of masses of the k-d point sets
	int k = 3; 
    // Size of point set 1 (transforming set) and point set 2 (fixed set);
    int nP1 = (int) samples.size();
    int nP2 = (int) mesh2->samples.size();

	float* m1 = new float[3];
    float* m2 = new float[3];
    // Calculate center of mesh1
	float sum[3] = {0.0f, 0.0f, 0.0f};
	for (int i = 0; i < nP1; i++)
		for (int c = 0; c < 3; c++)
			sum[c] += verts[ samples[i] ]->coords[c];
	for (int c = 0; c < 3; c++)
		m1[c] = sum[c] / (float) nP1;
    // Calculate center of mesh2
	sum[0] = sum[1] = sum[2] = 0.0f;
	for (int i = 0; i < nP2; i++)
		for (int c = 0; c < 3; c++)
			sum[c] += mesh2->verts[ mesh2->samples[i] ]->coords[c];
	for (int c = 0; c < 3; c++)
		m2[c] = sum[c] / (float) nP2;

	MatrixXf ms1(k, nP1), ms2(k, nP1); // Mean-shifted kxnP1 matrix for the mesh1.samples and matched samples in mesh2; X is for dynamic size f is for float
	Matrix3f cov12; // kxk cross-covariance mtrx of 2 sets, mesh & voxels where k = 3
	Matrix4f Q; //complicated matrix created using cov12

	float** R1 = new float*[3]; //desired 3x3 rotation mtrx of current iteration
	for (int i = 0; i < 3; i++)
		R1[i] = new float[3];
	float* t1 = new float[3]; //desired 1x3 translation

    // Parameters of the linear system to be solved in Adaptive Scaling
	float A, B;
    float* C = new float[3];
    float* D = new float[3];

	// Iteration parameters
	int nIters = 0;
	float mse = INF, prevMse = 0.0f; // Mean-square error objective function to be minimized
	float alignErr = 100.0f, closeness = 0.0000001f; //0.0000001f;//for face.xyz  //0.00001f;//for completehomer.xyz // 0.000000000001f default

    while (++nIters <= nMaxIters)
	{
        // break; //see the initial alignment (no ICP)
        // nIters > 70 && // Activate closeness test after a good deal of iters; a heuristic no in use
		if (fabs(mse - prevMse) <= closeness)
			break;
        
        // New center of mass of the transforming mesh
		sum[0] = sum[1] = sum[2] = 0.0f;
		for (int v = 0; v < nP1; v++)
			for (int c = 0; c < 3; c++)
				sum[c] += verts[ samples[v] ]->coords[c];
		for (int c = 0; c < 3; c++)
			m1[c] = sum[c] / (float) nP1;

        /*
            FIND CORRESPONDENCES
        */
        // true means 1-to-1 matching will be used between this mesh points and mesh2 points,
        // i.e., no mesh2 point appears more than once in the closest-pnt matching
        // (good for scaleAdaptive since it prevents selection of the same points over and over again when the scale difference is high)
		bool oneToOneMatches = oneToOne;
		int nMatched2 = 0; // Number of matched mesh2 vertices of the previous iteration

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

        float dist; //closest w.r.t. Euclidean distance
		vector< int > matchedP2s; //idxs of mesh2 points that are matched to mesh1 points

        //Brute force search for closest matches in quadratic time (n^2)
		for (int bv1 = 0; bv1 < nP1; bv1++) //for each mesh verts
		{
			float minDist = INF;
			int mv = -1;
			for (int bv2 = 0; bv2 < nP2; bv2++) //find its closest voxel match
			{
                dist = distanceBetween2(mesh2->verts[ mesh2->samples[bv2] ]->coords, verts[ samples[bv1] ]->coords); //omit sqrt() to save time
                if(dist > minDist)
                    continue;

				if (!oneToOneMatches) //previously selected point can still be selected
				{
                    minDist = dist;
                    mv = mesh2->samples[bv2];
				}
				else if (mesh2->verts[ mesh2->samples[bv2] ]->closestVertIdx == -1)	//previously selected point cannot be selected in this oneToOneMatches mode
				{
					minDist = dist;
                    mv = mesh2->samples[bv2];
				}
			}
			if (mv == -1)
				cout << "WARNING: this may happen in oneToOneMatches mode if a valid mesh2 point candidate cannot be found; easy fix: oneToOneMatches=false\n";
			matchedP2s.push_back(mv);
			mesh2->verts[mv]->closestVertIdx = samples[bv1];
		}

        // Now that closest pnt correspondence is done, i.e. closestVertIdx's are set, i can fill ms1
		// Fill mean-shifted kxnP1 ms1 matrix for the 1st set (mesh1)
		for (int i = 0; i < k; i++) // for all k rows
			for (int j = 0; j < nP1; j++) // fill columns
				ms1(i, j) = verts[ mesh2->verts[ matchedP2s[j] ]->closestVertIdx ]->coords[i] - m1[i];
        
        // Fill mean-shifted 3xnVerts matrix for the 2nd set (matched voxels, hence of size nVerts for sure)
		// New center of mass of the newly matched voxels
		sum[0] = sum[1] = sum[2] = 0.0f;
		for (int v = 0; v < nP1; v++)
			for (int c = 0; c < 3; c++)
				sum[c] += mesh2->verts[ matchedP2s[v] ]->coords[c];
		for (int c = 0; c < 3; c++)
			m2[c] = sum[c] / (float) nP1;
        
        //fill mean-shifted kxnP1 ms2 matrix for the mesh2
		for (int i = 0; i < k; i++) // For all k rows
			for (int j = 0; j < nP1; j++) // Fill n1 columns
				ms2(i, j) = mesh2->verts[ matchedP2s[j] ]->coords[i] - m2[i];
		MatrixXf ms2t = ms2.transpose(); //transpose of ms

        //kxk cross-covariance mtrx between mesh1 and mesh2
		cov12 = ms1 * ms2t; //k x nP1 * nP1 x k
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				cov12(i, j) /= nP1;

        //Q = 4x4 complicated mtrx, see the paper
		Q(0, 0) = cov12(0, 0) + cov12(1, 1) + cov12(2, 2);
		Q(1, 0) = Q(0, 1) = cov12(1, 2) - cov12(2, 1); //(paper notation: S23 - S32, say S to cov12, and idxs of paper, i.e. start from 1)
		Q(2, 0) = Q(0, 2) = cov12(2, 0) - cov12(0, 2); //S31 - S13
		Q(3, 0) = Q(0, 3) = cov12(0, 1) - cov12(1, 0); //S12 - S21
		Q(1, 1) = 2*cov12(0, 0) - Q(0, 0);		Q(1, 2) = cov12(0, 1) + cov12(1, 0);	Q(1, 3) = cov12(0, 2) + cov12(2, 0);
		Q(2, 1) = cov12(1, 0) + cov12(0, 1);	Q(2, 2) = 2*cov12(1, 1) - Q(0, 0);		Q(2, 3) = cov12(1, 2) + cov12(2, 1);
		Q(3, 1) = cov12(2, 0) + cov12(0, 2);	Q(3, 2) = cov12(2, 1) + cov12(1, 2);	Q(3, 3) = 2*cov12(2, 2) - Q(0, 0);

        // Desired rotation is decided by 1x4 quaternion vector, which is the (largest eigval's) eigvec of Q
		SelfAdjointEigenSolver<Matrix4f> eigensolver(Q); //solve eigenvalue/vector for symmetric n by n matrices, a.k.a. selfadjoint matrices
		if (eigensolver.info() != Success) abort();
		Vector4f eig_vals = eigensolver.eigenvalues(); // Given in ascending order
		Matrix4f eig_vecs = eigensolver.eigenvectors(); // Last column of eig_vecs corresponds to the largest eigenvalue, hence the one i want

        //1x4 quaternion = [q0, q1, q2, q3]
		float q0 = eig_vecs(0, 3), q1 = eig_vecs(1, 3), q2 = eig_vecs(2, 3), q3 = eig_vecs(3, 3);
		//if (q0 < 0) { q0 = -q0;	q1 = -q1;	q2 = -q2;	q3 = -q3; } //needed??????????

        // Fill desired 3x3 rotation mtrx (i'm used to 4x4 mtrxs w/ coords but this 3x3 works because of quaternions, i guess)
		R1[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;	R1[0][1] = 2*(q1*q2 - q0*q3);				R1[0][2] = 2*(q1*q3 + q0*q2);
		R1[1][0] = 2*(q1*q2 + q0*q3);				R1[1][1] = q0*q0 + q2*q2 - q1*q1 - q3*q3;	R1[1][2] = 2*(q2*q3 - q0*q1);
		R1[2][0] = 2*(q1*q3 - q0*q2);				R1[2][1] = 2*(q2*q3 + q0*q1);				R1[2][2] = q0*q0 + q3*q3 - q1*q1 - q2*q2;

        //our ICP that solves for translation and uniform scaling simulataneously to improve the initial scale guess as well

        //the math of this part is explained in my paper titled Skuller; check it out
        A = B = 0.0f;
        for (int i = 0; i < 3; i++)
            C[i] = D[i] = 0.0f;
        
        // Apply just ICP rotation of this iteration to all verts of the transforming mesh
        for (int v = 0; v < (int) verts.size(); v++)
        {
            float x = verts[v]->coords[0],
                    y = verts[v]->coords[1],
                    z = verts[v]->coords[2];
            //rotation
            verts[v]->coords[0] = R1[0][0]*x + R1[0][1]*y + R1[0][2]*z;
            verts[v]->coords[1] = R1[1][0]*x + R1[1][1]*y + R1[1][2]*z;
            verts[v]->coords[2] = R1[2][0]*x + R1[2][1]*y + R1[2][2]*z;
        }

        for (int v = 0; v < nP1; v++)
        {
            // D = [Dx Dy Dz] is computed Dx is the sum of x-coords of all matched mesh2 points,
            //   and so on (since matched mesh2 points change each iteration D is not fixed); 
            //   similarly C = [Cx Cy Cz] for transformed, i.e. ICP-rotated, mesh points
            for (int i = 0; i < 3; i++)
            {
                D[i] += mesh2->verts[ matchedP2s[v] ]->coords[i];
                C[i] += verts[ samples[v] ]->coords[i];
            }
            A += verts[ samples[v] ]->coords[0]*verts[ samples[v] ]->coords[0] + verts[ samples[v] ]->coords[1]*verts[ samples[v] ]->coords[1] + verts[ samples[v] ]->coords[2]*verts[ samples[v] ]->coords[2];
            B += verts[ samples[v] ]->coords[0]*mesh2->verts[ matchedP2s[v] ]->coords[0] + verts[ samples[v] ]->coords[1]*mesh2->verts[ matchedP2s[v] ]->coords[1] + verts[ samples[v] ]->coords[2]*mesh2->verts[ matchedP2s[v] ]->coords[2];
        }

        //solve for Px = r linear system to decide x = [s t1.x t1.y t1.z] where s is the uniform scale factor, t1[] is the translation
        Matrix4f P;
        P << A,C[0],C[1],C[2],	C[0],(float)nP1,0.0f,0.0f,	C[1],0.0f,(float)nP1,0.0f,	C[2],0.0f,0.0f,(float)nP1;
        //cout << "Here is the matrix P:\n" << P << endl;
        Vector4f r;
        r << B,D[0],D[1],D[2];
        //cout << "Here is the right hand side r:\n" << r << endl;
        Vector4f x = P.jacobiSvd(ComputeFullU | ComputeFullV).solve(r); //s = x(0) is the uniform scale factor, and the rest is the translation t1
        //cout << "The least-squares solution is:\n" << x << "\nsame as r?\n" << P * x << endl;

        //apply translation and scaling to all vertices
        for (int v = 0; v < (int) verts.size(); v++)
        {
            //scaling
            verts[v]->coords[0] *= x(0);
            verts[v]->coords[1] *= x(0);
            verts[v]->coords[2] *= x(0);

            //translation
            verts[v]->coords[0] += x(1);
            verts[v]->coords[1] += x(2);
            verts[v]->coords[2] += x(3);
        }
        //cout << x(0) << "\t\t" << x(1) << " " << x(2) << " " << x(3) << "\n";//cout << x(0) << " ";

        prevMse = mse; //mse of previous iteration as a termination condition
		//compute new mse value as an iteration condition
		mse = 0.0f;
		for (int bv2 = 0; bv2 < nP1; bv2++)
		{
			//squared distance b/w (target) matched voxel and its closest match in mesh whose coords just updated above
			dist = pow(mesh2->verts[ matchedP2s[bv2] ]->coords[0] - verts[ mesh2->verts[ matchedP2s[bv2] ]->closestVertIdx ]->coords[0], 2.0f) +
				   pow(mesh2->verts[ matchedP2s[bv2] ]->coords[1] - verts[ mesh2->verts[ matchedP2s[bv2] ]->closestVertIdx ]->coords[1], 2.0f) +
				   pow(mesh2->verts[ matchedP2s[bv2] ]->coords[2] - verts[ mesh2->verts[ matchedP2s[bv2] ]->closestVertIdx ]->coords[2], 2.0f);
			mse += dist; //no sqrt() above for efficiency, hence dist = squared dist indeed
		}
		mse /= nP1;
		if (nIters % 10 == 1)
			cout << "mse of ICP iter" << nIters-1 << ": " << mse << "\tmse diff: " << fabs(mse - prevMse) << endl;

    } // end of while nIters

    //recapture memo
	for (int i = 0; i < 3; i++) delete [] R1[i]; delete [] R1;
	delete [] t1;
	delete [] C;
	delete [] D;

    cout << nIters-1 << "'th ICP iteration w/ final mse: " << mse << "\nICP done!\n";

    return mse;
}

void Mesh::resultToFile(char* fName)
{
	//fprints transformed mesh
	FILE* fPtrS = fopen(fName, "w");
	for (int v = 0; v < (int) verts.size(); v++)
		fprintf(fPtrS, "%f %f %f\n", verts[v]->coords[0], verts[v]->coords[1], verts[v]->coords[2]);
	fclose(fPtrS);

	cout << "\ntransformed mesh fprinted to " << fName << endl;
}
