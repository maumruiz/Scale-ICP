#include <iostream>
#include "Mesh.h"
#include "icp.h"

using namespace std;

int main(int argc, char ** argv) {
    if (argc != 5)
	{
        cout << "4 arguments must have been provided (you did " << argc << "):\n" <<
			    "ScaleICP <landmarks-transforming> <landmarks-fixed> <mesh-transforming> <output-filepath>\n";
		exit(0);
    }

    int nMaxIters = 200;
	bool oneToOneStart = false;
    bool prescale = false;
    bool perfCorr = true;

	char fName[250], fName2[250], fName3[250], fNameExt[250], fNameExt2[250], fNameExt3[250];
    char fName4[250], ofName[250], ofNameExt[250], ofNameExt2[250];
	float minDisplacement = 10.0f;
    
    sprintf_s(fName, sizeof(fName), argv[1]);
    sprintf_s(fName2, sizeof(fName), argv[2]);
    sprintf_s(fName3, sizeof(fName), argv[3]);
    sprintf_s(ofName, sizeof(ofName), argv[4]);
	

    // // sprintf_s(fName, sizeof(fName), argv[1]);
    // // sprintf_s(fName2, sizeof(fName), argv[2]);
    // sprintf_s(fName, sizeof(fName), "luis_scaled_aligned_landmarks");
    // sprintf_s(fName2, sizeof(fName), "baseAligned_landmarks");
    // sprintf_s(fName3, sizeof(fName), "luis_scaled_aligned");

    sprintf_s(fNameExt, sizeof(fNameExt), "%s.xyz", fName);
    sprintf_s(fNameExt2, sizeof(fNameExt2), "%s.xyz", fName2);
    sprintf_s(fNameExt3, sizeof(fNameExt2), "%s.obj", fName3);

    Mesh* mesh1 = new Mesh(), * mesh2 = new Mesh(), * meshFinal = new Mesh();
	mesh1->loadPnt(fNameExt);
    mesh2->loadPnt(fNameExt2);
    meshFinal->loadObj(fNameExt3, false);

    // ICP* icp = new ICP();
    // icp->align(mesh1, mesh2, nMaxIters, oneToOneStart, minDisplacement, prescale);
    // mesh1->ICP(mesh2, nMaxIters, oneToOneStart, minDisplacement, prescale);
    vector<pair<Vector4f, float**>> transf = mesh1->ICP(mesh2, nMaxIters, oneToOneStart, minDisplacement, prescale, perfCorr);
    
	if (prescale) {
		for (int v = 0; v < (int)meshFinal->verts.size(); v++)
			for (int c = 0; c < 3; c++)
				meshFinal->verts[v]->coords[c] *= (mesh2->maxEucDist / mesh1->maxEucDist);
    }
    // icp->transform(meshFinal);

    for(int i=0; i<transf.size(); i++) {
        meshFinal->transform(transf[i].first, transf[i].second);
    }
    

    sprintf_s(ofNameExt, sizeof(ofNameExt), "%s_landmarks.xyz", ofName);
	mesh1->resultToFile(ofNameExt);

    sprintf_s(ofNameExt2, sizeof(ofNameExt), "%s.obj", ofName);
    meshFinal->updateObjVertices(fNameExt3, ofNameExt2);

    return 0;
}