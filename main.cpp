#include <iostream>
#include "Mesh.h"

using namespace std;

int main(int argc, char ** argv) {
    int nMaxIters = 200;
	bool oneToOneStart = true;

	char fName[250], fName2[250], fName3[250], fNameExt[250], fNameExt2[250], fNameExt3[250];
    char ofNameExt[250], ofNameExt2[250];
	float minDisplacement = 0.001f;

    // sprintf_s(fName, sizeof(fName), argv[1]);
    // sprintf_s(fName2, sizeof(fName), argv[2]);
    sprintf_s(fName, sizeof(fName), "alex_downsampled");
    sprintf_s(fName2, sizeof(fName), "targetHead");
    sprintf_s(fName3, sizeof(fName), "alex");

    sprintf_s(fNameExt, sizeof(fNameExt), "input\\%s.xyz", fName);
    sprintf_s(fNameExt2, sizeof(fNameExt2), "input\\%s.obj", fName2);
    sprintf_s(fNameExt3, sizeof(fNameExt2), "input\\%s.obj", fName3);

    Mesh* mesh1 = new Mesh(), * mesh2 = new Mesh(), * meshFinal = new Mesh();
	float mesh1Size = mesh1->loadPnt(fNameExt);
    mesh2->loadObj(fNameExt2, true);
    meshFinal->loadObj(fNameExt3, false);

    vector<pair<Vector4f, float**>> transf = mesh1->ICP(mesh2, nMaxIters, oneToOneStart, minDisplacement);

    sprintf_s(ofNameExt, sizeof(ofNameExt), "output\\%s_aligned.xyz", fName);
	mesh1->resultToFile(ofNameExt);

    sprintf_s(ofNameExt2, sizeof(ofNameExt), "output\\%s_aligned.obj", fName3);
    meshFinal->resultToObj(ofNameExt2);

    return 0;
}