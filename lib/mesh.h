#ifndef _MESH_H_
#define _MESH_H_

#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <vector>
#include <set>
using namespace std;

#include "Eigen/Dense"
using namespace Eigen;

#define INF (float) 1e38

struct Triangle
{
	int idx; //Index of the triangle
	int v1i, v2i, v3i; //Index of the vertices forming the triangle
	
	Triangle(int i, int i1, int i2, int i3) : idx(i), v1i(i1), v2i(i2), v3i(i3) {};
};

struct Vertex
{
	
	int idx; //Index of the vertex
	float* coords; //Coordinates
	int closestVertIdx; //Index of the closest vertex (for ICP)
	bool sample;

	Vertex(int i, float* c) : idx(i), coords(c), sample(false) {};
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
	vector< Triangle* > tris;
	vector< Vertex* > verts;
	vector< Edge* > edges;
	vector< int > samples;
	float maxEucDist;

	Mesh() {};
	
	float ICP(Mesh* mesh2, int nMaxIters, bool oneToOne, float minDisplacement);
	bool loadPnt(char* meshFile);
	void resultToFile(char* meshFile);

private:
	// int addVertex(float* v);
	// int addVertex(float x, float y, float z);
	void addVertex(float* c);
	// int addTriangle(int v1i, int v2i, int v3i);
	// void addEdge(int v1i, int v2i, bool interiorEdge = false);
};

#endif
