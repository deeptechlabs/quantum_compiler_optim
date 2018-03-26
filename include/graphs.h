#include "stdafx.h"


/* 
 * Data structures
 */

struct edge;
struct node;
struct node_cmp;
inline bool operator<(const edge& lhs, const edge& rhs);


/*
 * Global variables
 */

int** dist;
int positions;

unsigned long ngates;
unsigned int nqubits;

std::set<edge> graph;
vector<vector<gate> > layers;
priority_queue<node, std::vector<node>, node_cmp> nodes;


/*
 * Functions here
 */

void build_graph_QX3();
bool contains(vector<int> v, int e);
int bfs(int start, int goal, std::set<edge>& graph);
void build_dist_table(std::set<edge>& graph);
