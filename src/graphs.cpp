#include "stdafx.h"
#include "graphs.h"

using namespace std;

namespace graphs {

	/*
	 * Global variables
	 */

	int** dist;
	int positions;

	unsigned long ngates = 0;
	unsigned int nqubits = 0;

	std::set<edge> graph;
	vector<vector<gate> > layers;
	priority_queue<node, std::vector<node>, node_cmp> nodes;

	//build a graph representing the coupling map of IBM QX3
	void build_graph_QX3() {
		graph.clear();
		positions = 16;
		edge e;
		e.v1 = 0;
		e.v2 = 1;
		graph.insert(e);
		e.v1 = 1;
		e.v2 = 2;
		graph.insert(e);
		e.v1 = 2;
		e.v2 = 3;
		graph.insert(e);
		e.v1 = 3;
		e.v2 = 14;
		graph.insert(e);
		e.v1 = 4;
		e.v2 = 3;
		graph.insert(e);
		e.v1 = 4;
		e.v2 = 5;
		graph.insert(e);
		e.v1 = 6;
		e.v2 = 7;
		graph.insert(e);
		e.v1 = 6;
		e.v2 = 11;
		graph.insert(e);
		e.v1 = 7;
		e.v2 = 10;
		graph.insert(e);
		e.v1 = 8;
		e.v2 = 7;
		graph.insert(e);
		e.v1 = 9;
		e.v2 = 8;
		graph.insert(e);
		e.v1 = 9;
		e.v2 = 10;
		graph.insert(e);
		e.v1 = 11;
		e.v2 = 10;
		graph.insert(e);
		e.v1 = 12;
		e.v2 = 5;
		graph.insert(e);
		e.v1 = 12;
		e.v2 = 11;
		graph.insert(e);
		e.v1 = 12;
		e.v2 = 13;
		graph.insert(e);
		e.v1 = 13;
		e.v2 = 4;
		graph.insert(e);
		e.v1 = 13;
		e.v2 = 14;
		graph.insert(e);
		e.v1 = 15;
		e.v2 = 0;
		graph.insert(e);
		e.v1 = 15;
		e.v2 = 14;
		graph.insert(e);
	}

	bool contains(vector<int> v, int e) {
		for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
			if (*it == e) {
				return true;
			}
		}
		return false;
	}

	//Breadth first search algorithm to determine the shortest paths between two physical qubits
	int bfs(int start, int goal, set<edge>& graph) {
		queue<vector<int> > queue;
		vector<int> v;
		v.push_back(start);
		queue.push(v);
		vector<vector<int> > solutions;

		int length;
		std::set<int> successors;
		while (!queue.empty()) {
			v = queue.front();
			queue.pop();
			int current = v[v.size() - 1];
			if (current == goal) {
				length = v.size();
				solutions.push_back(v);
				break;
			} else {
				successors.clear();
				for (set<edge>::iterator it = graph.begin(); it != graph.end();
						it++) {
					edge e = *it;
					if (e.v1 == current && !contains(v, e.v2)) {
						successors.insert(e.v2);
					}
					if (e.v2 == current && !contains(v, e.v1)) {
						successors.insert(e.v1);
					}
				}
				for (set<int>::iterator it = successors.begin();
						it != successors.end(); it++) {
					vector<int> v2 = v;
					v2.push_back(*it);
					queue.push(v2);
				}
			}
		}
		while (!queue.empty() && queue.front().size() == length) {
			if (queue.front()[queue.front().size() - 1] == goal) {
				solutions.push_back(queue.front());
			}
			queue.pop();
		}

		for (int i = 0; i < solutions.size(); i++) {
			vector<int> v = solutions[i];
			for (int j = 0; j < v.size() - 1; j++) {
				edge e;
				e.v1 = v[j];
				e.v2 = v[j + 1];
				if (graph.find(e) != graph.end()) {
					return (length-2)*7;
					return length - 2;
				}
			}
		}

		return (length - 2)*7 + 4;
	}

	void build_dist_table(std::set<edge>& graph) {
		dist = new int*[positions];

		for (int i = 0; i < positions; i++) {
			dist[i] = new int[positions];
		}

		for (int i = 0; i < positions; i++) {
			for (int j = 0; j < positions; j++) {
				if (i != j) {
					dist[i][j] = bfs(i,j,graph);
				} else {
					dist[i][i] = 0;
				}
			}
		}
	}

}
