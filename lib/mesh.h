#ifndef _MESH_H_
#define _MESH_H_

#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <fstream>
using namespace std;

#include "Eigen/Dense"
using namespace Eigen;

#define INF (float) 1e38

struct Face
{
	int idx; //Index of the triangle
	int* vs; //Index of the vertices forming the triangle
	
	Face(int i, int* vis) : idx(i), vs(vis) {};
};

struct Vertex
{
	
	int idx; //Index of the vertex
	float* coords; //Coordinates
	int closestVertIdx; //Index of the closest vertex (for ICP)

	Vertex(int i, float* c) : idx(i), coords(c) {};
};

struct Edge
{
	int idx; //Index of the edge
	int v1i, v2i; //Index of vertices forming the edge
	float length; //Distance between v1i & v2i
	bool interior; //true if this is an interior edge, false if it is on the boundary (i.e. useful for geodesic computations)

	Edge(int i, int i1, int i2, float l, bool in) : idx(i), v1i(i1), v2i(i2), length(l), interior(in) {};
};

class Mesh
{
public:
	vector< Face* > faces;
	vector< Vertex* > verts;
	vector< Edge* > edges;
	vector< int > samples;
	int numVerts;
	int numFaces;
	float maxEucDist;

	Mesh() {};
	
	vector<pair<Vector4f, float**>> ICP(Mesh* mesh2, int nMaxIters, bool oneToOne, float minDisplacement, bool prescale, bool perfCorr);
	void transform(Vector4f, float**);
	bool loadPnt(char* meshFile);
	bool loadObj(char* meshFile, bool calculateSize);
	bool loadLandmarks(char* meshFile);
	void resultToFile(char* meshFile);
	void resultToObj(char* meshFile);
	void updateObjVertices(char* meshFile, char* outFile);

private:
	// int addVertex(float* v);
	// int addVertex(float x, float y, float z);
	float calculateMaxEucDist();
	void addVertex(float* c);
	int addFace(int* fa);
	// void addEdge(int v1i, int v2i, bool interiorEdge = false);
};

#endif
