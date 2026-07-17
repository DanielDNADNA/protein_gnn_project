#pragma once
#include <vector>

struct Edges {

	std::vector<int> sourceNode; 
	std::vector<int> sinkNode;
	std::vector<int> edgeType;
	std::vector<float> euclideanDistance;

	void reserve(size_t n) {

		sourceNode.reserve(6 * n);
		sinkNode.reserve(6 * n);
		edgeType.reserve(6 * n);
		euclideanDistance.reserve(6 * n);
	}
};