#include <iostream>
#include "Mesh.h"

using namespace std;

int main(int argc, char ** argv) {
    int nMaxIters = 200;
	bool oneToOneStart = true;

	char fName[250], fName2[250], fNameExt[250], fNameExt2[250], ofNameExt[250], ofNameExt2[250];
	float minDisplacement = 0.001f;

    // sprintf_s(fName, sizeof(fName), "targetHead");
    // sprintf_s(fNameExt, sizeof(fNameExt), "input\\%s.obj", fName);
    // Mesh* mesh1 = new Mesh();
    // mesh1->loadObj(fNameExt);
    // sprintf_s(ofNameExt, sizeof(ofNameExt), "output\\%s_test.obj", fName);
    // mesh1->resultToObj(ofNameExt);

    // sprintf_s(fName, sizeof(fName), argv[1]);
    // sprintf_s(fName2, sizeof(fName), argv[2]);
    sprintf_s(fName, sizeof(fName), "st");
    sprintf_s(fName2, sizeof(fName), "targetHead");

    sprintf_s(fNameExt, sizeof(fNameExt), "input\\%s.obj", fName);
    sprintf_s(fNameExt2, sizeof(fNameExt2), "input\\%s.obj", fName2);

    Mesh* mesh1 = new Mesh(), * mesh2 = new Mesh(), * meshFinal = new Mesh();
	mesh1->loadObj(fNameExt);
    mesh2->loadObj(fNameExt2);
    // meshFinal->loadObj(fNameExt);

    vector<pair<Vector4f, float**>> transf = mesh1->ICP(mesh2, nMaxIters, oneToOneStart, minDisplacement);


    // for(int i=0; i<transf.size(); i++) {
    //     meshFinal->transform(transf[i].first, transf[i].second);
    // }

    sprintf_s(ofNameExt, sizeof(ofNameExt), "output\\%s_aligned.obj", fName);
	mesh1->resultToObj(ofNameExt);

    // sprintf_s(ofNameExt2, sizeof(ofNameExt), "output\\%s_aligned_2.obj", fName);
    // meshFinal->resultToObj(ofNameExt2);

    return 0;
}