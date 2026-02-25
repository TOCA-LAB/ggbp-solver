#ifndef DECODER_H
#define DECODER_H

#include <vector>
#include <forward_list>
#include <algorithm>
#include "ggbp/graph.hpp"

using namespace std;

class Decoder {
public:

	Decoder(Graph &g, int const _r);

	~Decoder();

	int countBurningNeighbours(const int &v, const vector<int> &burned) const;

	tuple<int, vector<int>> decoder(const vector< double >& chromosome) const;

	double decode(const vector< double >& chromosome) const;

	vector<int> getBurningSequence(const std::vector< double >& chromosome) const;

private:

	int r;
    Graph graph;
	int n_vertices;
};

#endif
