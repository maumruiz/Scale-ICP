#include <iostream>
#include "Mesh.h"

using namespace std;

int main(int argc, char ** argv) {
    int nMaxIters = 200;
	bool oneToOneStart = true;

	char fName[250], fName2[250], fNameExt[250], fNameExt2[250], ofNameExt[250];
	float minDisplacement = 0.001f;


    // sprintf_s(fName, sizeof(fName), argv[1]);
    // sprintf_s(fName2, sizeof(fName), argv[2]);
    sprintf_s(fName, sizeof(fName), "rt1");
    sprintf_s(fName2, sizeof(fName), "baseMesh");

    sprintf_s(fNameExt, sizeof(fNameExt), "input\\%s.xyz", fName);
    sprintf_s(fNameExt2, sizeof(fNameExt2), "input\\%s.xyz", fName2);

    Mesh* mesh1 = new Mesh(), * mesh2 = new Mesh();
	mesh1->loadPnt(fNameExt);
    mesh2->loadPnt(fNameExt2);

    mesh1->ICP(mesh2, nMaxIters, oneToOneStart, minDisplacement);

    sprintf_s(ofNameExt, sizeof(ofNameExt), "output\\%s_aligned.xyz", fName);
	mesh1->resultToFile(ofNameExt);

    return 0;
}