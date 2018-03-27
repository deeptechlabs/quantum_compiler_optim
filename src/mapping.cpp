#include "mapping.h"

using namespace std;

//A very simplified QASM parser
void read_qasm(std::ifstream& infile) {
	std::string line;

	std::getline(infile, line);
	if(line != "OPENQASM 2.0;") {
		cerr << "ERROR: first line of the file has to be: OPENQASM 2.0;" << endl;
		exit(1);
	}

	std::getline(infile, line);
	if(line != "include \"qelib1.inc\";") {
		cerr << "ERROR: second line of the file has to be: include \"qelib1.inc\"" << endl;
		exit(1);
	}

	std::getline(infile, line);
	int n = -1;
	if (sscanf(line.c_str(), "qreg q[%d];", &n) != 1) {
		cerr << "ERROR: failed to parse qasm file: " << line << endl;
		exit(1);
	}
	if (n > positions) {
		cerr << "ERROR: too many qubits for target architecture: " << n << endl;
		exit(2);
	}

	std::getline(infile, line);

	int* last_layer = new int[n];
	for (int i = 0; i < n; i++) {
		last_layer[i] = -1;
	}

	while (std::getline(infile, line)) {

		if (line == "") {
			continue;
		}
		gate g;
		int layer;

		int nq = sscanf(line.c_str(), "%3s q[%d],q[%d];", g.type, &g.control,
				&g.target);

		if (nq == 3) {
			layer = max(last_layer[g.target], last_layer[g.control]) + 1;
			last_layer[g.target] = last_layer[g.control] = layer;
		} else if (nq == 2) {
			g.target = g.control;
			g.control = -1;
			layer = last_layer[g.target] + 1;
			last_layer[g.target] = layer;
		} else {
			double angle;
			if(sscanf(line.c_str(), "rz(%f) q[%d];", &angle, &g.target) == 2) {
				g.control = -1;
				strcpy(g.type, "rz");
				layer = last_layer[g.target] + 1;
				last_layer[g.target] = layer;
			} else {
				cerr << "ERROR: could not read gate: " << line << endl;
				exit(1);
			}
		}
		ngates++;

		if (layers.size() <= layer) {
			layers.push_back(vector<gate>());
		}
		layers[layer].push_back(g);
	}

	nqubits = n;
	delete[] last_layer;
}

void expand_node(const vector<int>& qubits, int qubit, edge *swaps, int nswaps,
		int* used, node base_node, const vector<gate>& gates, int** dist, int next_layer) {

	if (qubit == qubits.size()) {
		//base case: insert node into queue
		if (nswaps == 0) {
			return;
		}
		node new_node;

		new_node.qubits = new int[positions];
		new_node.locations = new int[nqubits];

		memcpy(new_node.qubits, base_node.qubits, sizeof(int) * positions);
		memcpy(new_node.locations, base_node.locations, sizeof(int) * nqubits);

		new_node.swaps = vector<vector<edge> >();
		new_node.nswaps = base_node.nswaps + nswaps;
		for (vector<vector<edge> >::iterator it2 = base_node.swaps.begin();
				it2 != base_node.swaps.end(); it2++) {
			vector<edge> new_v(*it2);
			new_node.swaps.push_back(new_v);
		}

		new_node.depth = base_node.depth + 5;
		new_node.cost_fixed = base_node.cost_fixed + 7 * nswaps;
		new_node.cost_heur = 0;

		vector<edge> new_swaps;
		for (int i = 0; i < nswaps; i++) {
			new_swaps.push_back(swaps[i]);
			int tmp_qubit1 = new_node.qubits[swaps[i].v1];
			int tmp_qubit2 = new_node.qubits[swaps[i].v2];

			new_node.qubits[swaps[i].v1] = tmp_qubit2;
			new_node.qubits[swaps[i].v2] = tmp_qubit1;

			if (tmp_qubit1 != -1) {
				new_node.locations[tmp_qubit1] = swaps[i].v2;
			}
			if (tmp_qubit2 != -1) {
				new_node.locations[tmp_qubit2] = swaps[i].v1;
			}
		}
		new_node.swaps.push_back(new_swaps);
		new_node.done = 1;

		for (vector<gate>::const_iterator it = gates.begin(); it != gates.end();
				it++) {
			gate g = *it;
			if (g.control != -1) {
				new_node.cost_heur = new_node.cost_heur + dist[new_node.locations[g.control]][new_node.locations[g.target]];
				if(dist[new_node.locations[g.control]][new_node.locations[g.target]] > 4) {
					new_node.done = 0;
				}
			}
		}

		//Calculate heuristics for the cost of the following layer
		new_node.cost_heur2 = 0;
		if(next_layer != -1) {
			for (vector<gate>::const_iterator it = layers[next_layer].begin(); it != layers[next_layer].end();
							it++) {
				gate g = *it;
				if (g.control != -1) {
					if(new_node.locations[g.control] == -1 && new_node.locations[g.target]) {
					//ignore this case
					} else if(new_node.locations[g.control] == -1) {
						int min = 1000;
						for(int i=0; i< positions; i++) {
							if(new_node.qubits[i] == -1 && dist[i][new_node.locations[g.target]] < min) {
								min = dist[i][new_node.locations[g.target]];
							}
						}
						new_node.cost_heur2 = new_node.cost_heur2 + min;
					} else if(new_node.locations[g.target] == -1) {
						int min = 1000;
						for(int i=0; i< positions; i++) {
							if(new_node.qubits[i] == -1 && dist[new_node.locations[g.control]][i] < min) {
								min = dist[new_node.locations[g.control]][i];
							}
						}
						new_node.cost_heur2 = new_node.cost_heur2 + min;
					} else {
						new_node.cost_heur2 = new_node.cost_heur2 + dist[new_node.locations[g.control]][new_node.locations[g.target]];
					}
				}
			}
		}

		nodes.push(new_node);
	} else {
		expand_node(qubits, qubit + 1, swaps, nswaps, used, base_node, gates,
				dist, next_layer);

		for (set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
			edge e = *it;
			if (e.v1 == base_node.locations[qubits[qubit]]
					|| e.v2 == base_node.locations[qubits[qubit]]) {
				if (!used[e.v1] && !used[e.v2]) {
					used[e.v1] = 1;
					used[e.v2] = 1;
					swaps[nswaps].v1 = e.v1;
					swaps[nswaps].v2 = e.v2;
					expand_node(qubits, qubit + 1, swaps, nswaps + 1, used,
							base_node, gates, dist, next_layer);
					used[e.v1] = 0;
					used[e.v2] = 0;
				}
			}
		}
	}
}

int getNextLayer(int layer) {
	int next_layer = layer+1;
	while(next_layer < layers.size()) {
		for(vector<gate>::iterator it = layers[next_layer].begin(); it != layers[next_layer].end(); it++) {
			if(it->control != -1) {
				return next_layer;
			}
		}
		next_layer++;
	}
	return -1;
}

node a_star_fixlayer(int layer, int* map, int* loc, int** dist) {

	int next_layer = getNextLayer(layer);

	node n;
	n.cost_fixed = 0;
	n.cost_heur = n.cost_heur2 = 0;
	n.qubits = new int[positions];
	n.locations = new int[nqubits];
	n.swaps = vector<vector<edge> >();
	n.done = 1;

	vector<gate> v = vector<gate>(layers[layer]);
	vector<int> considered_qubits;

	//Find a mapping for all logical qubits in the CNOTs of the layer that are not yet mapped
	for (vector<gate>::iterator it = v.begin(); it != v.end(); it++) {
		gate g = *it;
		if (g.control != -1) {
			considered_qubits.push_back(g.control);
			considered_qubits.push_back(g.target);
			if(loc[g.control] == -1 && loc[g.target] == -1) {
				set<edge> possible_edges;
				for(set<edge>::iterator it = graph.begin(); it != graph.end(); it++) {
					if(map[it->v1] == -1 && map[it->v2] == -1) {
						possible_edges.insert(*it);
					}
				}
				if(!possible_edges.empty()) {
					edge e = *possible_edges.begin();
					loc[g.control] = e.v1;
					map[e.v1] = g.control;
					loc[g.target] = e.v2;
					map[e.v2] = g.target;
				} else {
					cout << "no edge available!";
					exit(1);
				}
			} else if(loc[g.control] == -1) {
				int min = 1000;
				int min_pos = -1;
				for(int i=0; i< positions; i++) {
					if(map[i] == -1 && dist[i][loc[g.target]] < min) {
						min = dist[i][loc[g.target]];
						min_pos = i;
					}
				}
				map[min_pos] = g.control;
				loc[g.control] = min_pos;
			} else if(loc[g.target] == -1) {
				int min = 1000;
				int min_pos = -1;
				for(int i=0; i< positions; i++) {
					if(map[i] == -1 && dist[loc[g.control]][i] < min) {
						min = dist[loc[g.control]][i];
						min_pos = i;
					}
				}
				map[min_pos] = g.target;
				loc[g.target] = min_pos;
			}
			n.cost_heur = max(n.cost_heur, dist[loc[g.control]][loc[g.target]]);
		} else {
		}
	}

	if(n.cost_heur > 4) {
		n.done = 0;
	}

	memcpy(n.qubits, map, sizeof(int) * positions);
	memcpy(n.locations, loc, sizeof(int) * nqubits);

	nodes.push(n);

	int *used = new int[positions];
	for (int i = 0; i < positions; i++) {
		used[i] = 0;
	}
	edge *edges = new edge[considered_qubits.size()];

	//Perform an A* search to find the cheapest permuation
	while (!nodes.top().done) {
		node n = nodes.top();
		nodes.pop();

		expand_node(considered_qubits, 0, edges, 0, used, n, v, dist, next_layer);

		delete[] n.locations;
		delete[] n.qubits;
	}

	node result = nodes.top();
	nodes.pop();

	//clean up
	delete[] used;
	delete[] edges;
	while (!nodes.empty()) {
		node n = nodes.top();
		nodes.pop();
		delete[] n.locations;
		delete[] n.qubits;
	}
	return result;
}
