#include "stdafx.h"
#include "graphs.h"

using namespace std;

namespace mapping {

	void read_qasm(std::ifstream& infile);

	void expand_node(const vector<int>& qubits, int qubit, edge *swaps, int nswaps, int* used, node base_node, const vector<gate>& gates, int** dist, int next_layer);


	int getNextLayer(int layer);

	node a_star_fixlayer(int layer, int* map, int* loc, int** dist);

}
