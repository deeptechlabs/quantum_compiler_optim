#include "stdafx.h"
#include "rl.h"
#include "graphs.h"
#include "mapping.h"

using namespace mapping;
using namespace graphs;

int main(int argc, char** argv) {

	if(argc != 3) {
		cout << "Usage: ./imb_mapping <input_file> <output_file>" << endl;
		exit(0);
	}

	graphs::build_graph_QX3();
	graphs::build_dist_table(graph);

	ifstream infile(argv[1]);
	mapping::read_qasm(infile);

	unsigned int width = 0;
	for (vector<vector<gate> >::iterator it = layers.begin(); it != layers.end(); it++) {
		if ((*it).size() > width) {
			width = (*it).size();
		}
	}

	char* bName = argv[1];
	cout << "Circuit name: " << bName << " (requires " << nqubits << " qubits)" << endl;

	cout << endl << "Before mapping: " << endl;
	cout << "  elementary gates: " << ngates << endl;
	cout << "  depth: " << layers.size() << endl;

	int *locations = new int[nqubits];
	int *qubits = new int[positions];

	for (int i = 0; i < positions; i++) {
		qubits[i] = -1;
	}
	for(int i = 0; i < nqubits; i++) {
		locations[i] = qubits[i] = i;
	}

	//Start mapping algorithm
	clock_t begin_time = clock();

	//Initially, no physical qubit is occupied
	for (int i = 0; i < positions; i++) {
			qubits[i] = -1;
	}

	//Initially, no logical qubit is mapped to a physical one
	for(int i = 0; i < nqubits; i++) {
		locations[i] = -1;
	}

	vector<gate> all_gates;
	int total_swaps = 0;

	//Fix the mapping of each layer
	for (int i = 0; i < layers.size(); i++) {
		node result = mapping::a_star_fixlayer(i, qubits, locations, dist);

		delete[] locations;
		delete[] qubits;
		locations = result.locations;
		qubits = result.qubits;

		vector<gate> h_gates = vector<gate>();

		//The first layer does not require a permutation of the qubits
		if (i != 0) {
			//Add the required SWAPs to the circuits
			for (vector<vector<edge> >::iterator it = result.swaps.begin();
					it != result.swaps.end(); it++) {
				for (vector<edge>::iterator it2 = it->begin(); it2 != it->end();
						it2++) {

					edge e = *it2;
					gate cnot;
					gate h1;
					gate h2;
					if (graph.find(e) != graph.end()) {
						cnot.control = e.v1;
						cnot.target = e.v2;
					} else {
						cnot.control = e.v2;
						cnot.target = e.v1;

						int tmp = e.v1;
						e.v1 = e.v2;
						e.v2 = tmp;
						if (graph.find(e) == graph.end()) {
							cerr << "ERROR: invalid SWAP gate" << endl;
							exit(2);
						}
					}
					strcpy(cnot.type, "cx");
					strcpy(h1.type, "h");
					strcpy(h2.type, "h");
					h1.control = h2.control = -1;
					h1.target = e.v1;
					h2.target = e.v2;

					gate gg;
					gg.control = cnot.control;
					gg.target = cnot.target;
					strcpy(gg.type, "SWP");

					all_gates.push_back(cnot);
					all_gates.push_back(h1);
					all_gates.push_back(h2);
					all_gates.push_back(cnot);
					all_gates.push_back(h1);
					all_gates.push_back(h2);
					all_gates.push_back(cnot);
					//Insert a dummy SWAP gate to allow for tracking the positions of the logical qubits
					all_gates.push_back(gg);
					total_swaps++;
				}
			}
		}

		//Add all gates of the layer to the circuit
		vector<gate> layer_vec = layers[i];
		for (vector<gate>::iterator it = layer_vec.begin();
				it != layer_vec.end(); it++) {
			gate g = *it;
			if (g.control == -1) {
				//single qubit gate
				if(locations[g.target] == -1) {
					//handle the case that the qubit is not yet mapped. This happens if the qubit has not yet occurred in a CNOT gate
					gate g2 = g;
					g2.target = -g.target -1;
					all_gates.push_back(g2);
				} else {
					//Add the gate to the circuit
					g.target = locations[g.target];
					all_gates.push_back(g);
				}
			} else {
				//CNOT gate
				g.target = locations[g.target];
				g.control = locations[g.control];

				edge e;
				e.v1 = g.control;
				e.v2 = g.target;

				if (graph.find(e) == graph.end()) {
					//flip the direction of the CNOT by inserting H gates
					e.v1 = g.target;
					e.v2 = g.control;
					if (graph.find(e) == graph.end()) {
						cerr << "ERROR: invalid CNOT: " << e.v1 << " - " << e.v2
								<< endl;
						exit(3);
					}
					gate h;
					h.control = -1;
					strcpy(h.type, "h");
					h.target = g.target;
					all_gates.push_back(h);

					h_gates.push_back(h);
					h.target = g.control;
					all_gates.push_back(h);

					h_gates.push_back(h);
					int tmp = g.target;
					g.target = g.control;
					g.control = tmp;
				}
				all_gates.push_back(g);
			}
		}
		if (h_gates.size() != 0) {
			if (result.cost_heur == 0) {
				cerr << "ERROR: invalid heuristic cost!" << endl;
				exit(2);
			}

			for (vector<gate>::iterator it = h_gates.begin();
					it != h_gates.end(); it++) {
				all_gates.push_back(*it);
			}
		}

	}

	//Fix the position of the single qubit gates
	for(vector<gate>::reverse_iterator it = all_gates.rbegin(); it != all_gates.rend(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];
			qubits[it->control] = tmp_qubit2;
			qubits[it->target] = tmp_qubit1;

			if(tmp_qubit1 != -1) {
				locations[tmp_qubit1] = it->target;
			}
			if(tmp_qubit2 != -1) {
				locations[tmp_qubit2] = it->control;
			}
		}
		if(it->target < 0) {
			int target = -(it->target +1);
			it->target = locations[target];
			if(locations[target] == -1) {
				//This qubit occurs only in single qubit gates -> it can be mapped to an arbirary physical qubit
				int loc = 0;
				while(qubits[loc] != -1) {
					loc++;
				}
				locations[target] = loc;
			}
		}
	}


	int *last_layer = new int[positions];
	for(int i=0; i<positions; i++) {
		last_layer[i] = -1;
	}

	vector<vector<gate> > mapped_circuit;


	//build resulting circuit
	for(vector<gate>::iterator it = all_gates.begin(); it != all_gates.end(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			continue;
		}
		if(it->control == -1) {
			//single qubit gate
			gate g = *it;
			int layer = last_layer[g.target] + 1;

			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<gate>());
			}
			mapped_circuit[layer].push_back(g);
			last_layer[g.target] = layer;
		} else {
			gate g = *it;
			int layer = max(last_layer[g.control], last_layer[g.target]) + 1;
			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<gate>());
			}
			mapped_circuit[layer].push_back(g);

			last_layer[g.target] = layer;
			last_layer[g.control] = layer;
		}
	}

	double time = double(clock() - begin_time) / CLOCKS_PER_SEC;

	cout << endl << "After mapping (no post mapping optimizations are conducted): " << endl;
	cout << "  elementary gates: " << all_gates.size()-total_swaps << endl;
	cout << "  depth: " << mapped_circuit.size() << endl;

	cout << endl << "The mapping required " << time << " seconds" << endl;

	cout << endl << "Initial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX3 architecture: " << endl;

	for(int i=0; i<nqubits; i++) {
		cout << "  q" << i << " is initially mapped to Q" << locations[i] << endl;
	}

	//Dump resulting circuit

	ofstream of(argv[2]);

	of << "OPENQASM 2.0;" << endl;
	of << "include \"qelib1.inc\";" << endl;
	of << "qreg q[16];" << endl;
	of << "creg c[16];" << endl;

	for (vector<vector<gate> >::iterator it = mapped_circuit.begin();
			it != mapped_circuit.end(); it++) {
		vector<gate> v = *it;
		for (vector<gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
			of << it2->type << " ";
			if (it2->control != -1) {
				of << "q[" << it2->control << "],";
			}
			of << "q[" << it2->target << "];" << endl;
		}
	}

	delete[] locations;
	delete[] qubits;
	delete[] last_layer;

	return 0;
}
