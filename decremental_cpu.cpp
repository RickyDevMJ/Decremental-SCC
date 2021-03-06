#include <queue>
#include <unordered_map>
#include <unordered_set>

// data-structure to store dynamic graph
class dyn_graph {

	public:

		// number of vertices
		long n;
		// csr representation of in-edges and out-edges
		long *in_row; 
		long *out_row;
		long *in_col; 
		long *out_col;
		// in degree and out degree for each vertex
		long *in_deg;
		long *out_deg;

		// maps vertex in current node to vertex in the base graph
		std::unordered_map<long, long> vertex_map;

		// function to remove an edge from the graph
		void remove_edge(long src, long dst) {
			#pragma omp parallel for
			for (long i = out_row[src]; i < out_row[src + 1]; i++) {
				if (out_col[i] == dst) {
					out_col[i] = -1;

					#pragma omp atomic
					out_deg[src]--;
				}
			}

			#pragma omp parallel for
			for (long i = in_row[dst]; i < in_row[dst + 1]; i++) {
				if (in_col[i] == src) {
					in_col[i] = -1;

					#pragma omp atomic
					in_deg[dst]--;
				}
			}
		}

		// function to remove a node from the graph
		void remove_node(long v) {
			// todo
		}
};

// data-structure to store a node of the SCC-Tree
class scc_tree {
	public:

		// unique identifier for every node
		long id;
		// The graph represented by the node, NULL if it is leaf
		dyn_graph G;
		// parent pointer of the node in the tree
		scc_tree *parent;
		// vector like representation of children of current node
		long max_children, num_children;
		scc_tree **children;
		long size;
		// this array stores all the vertices of base graph contained in the current node
		long *vertices;
};

// simple structure to hold the unreachable vertices and edges; return value of find_unreachable
struct unreachable
{
	std::vector<long> U, I_src, I_dst;
	std::unordered_set<long> U_set;
};

// simple function to calaculate depth of a node from root of the tree
long depth(scc_tree *node)
{
	long d = -1;
	while (node)
	{
		d++;
		node = node->parent;
	}
	return d;
}

// simple function to find least common ancestor of two nodes in the tree
scc_tree *LCA(scc_tree *n1, scc_tree *n2)
{
	long d1 = depth(n1), d2 = depth(n2);
	long diff = d1 - d2;

	// If n2 is deeper, swap n1 and n2
	if (diff < 0)
	{
		scc_tree *temp = n1;
		n1 = n2;
		n2 = temp;
		diff = -diff;
	}

	// Move n1 up until it reaches the same level as n2
	while (diff--)
		n1 = n1->parent;

	while (n1 && n2)
	{
		if (n1 == n2)
			return n1;
		n1 = n1->parent;
		n2 = n2->parent;
	}

	return NULL;
}

// helper function to recursively lift up the deleted node in the SCC-tree
void lift_up(scc_tree *T, unreachable &R, scc_tree **t_array) {

	// if node is the root
	if (T->parent == NULL) {
		#pragma omp parallel for
		for (long i = 0; i < T->max_children; i++) {
			if (T->children[i] == NULL)
				continue;

			for (long j = 0; j < R.U.size(); j++) {
				if (T->children[i] == t_array[R.U[j]]) {
					t_array[R.U[j]]->parent = NULL;

					T->children[i] = NULL;
					#pragma omp atomic
					T->num_children--;

					// todo : modify T.G

					break;
				}
			}
		}
	}
	else {
		#pragma omp parallel for
		for (long i = 0; i < T->max_children; i++) {
			if (T->children[i] == NULL)
				continue;

			for (long j = 0; j < R.U.size(); j++) {
				if (T->children[i] == t_array[R.U[j]]) {
					t_array[R.U[j]]->parent = T->parent;

					T->children[i] = NULL;
					#pragma omp atomic
					T->num_children--;

					// todo : modify T.G and T->parent and T->parent.G

					break;
				}
			}
		}

		// todo: recursive call
	}
}

// function to remove an edge from tree and subsequently from the base graph
void remove_edge(long src, long dst, long n, scc_tree **t_array) {
	scc_tree *T = LCA(t_array[src], t_array[dst]);

	// edge across two different SCCs
	if (T == NULL) {
		// todo
		return;
	}

	long S[2];
	S[0] = T->G.vertex_map[src];
	S[1] = T->G.vertex_map[dst];
	unreachable R = find_unreachable_down(T->G, 2, S, T->G.vertex_map[T->id]);
	unreachable R1 = find_unreachable_up(T->G, 2, S, T->G.vertex_map[T->id]);

	for (long i = 0; i < R1.U.size(); i++) {
		if (R.U_set.find(R1.U[i]) == R.U_set.end()) {
			R.U.push_back(R1.U[i]);
			// todo : update R.I_src and R.I_dst
			R.U_set.insert(R1.U[i]);	
		}
	}

	lift_up(T, R, t_array);

}

// function to find the set of vertices not reachable from the source w
unreachable find_unreachable_down(dyn_graph G, long ns, long *S, long w) {
	std::queue<long> Q[2];

	for (long i = 0; i < ns; i++) {
		long v = S[i];
		if (v != w) {
			// if v has no in-edges
			if (G.in_deg[v] == 0) {
				Q[0].push(v);
			}
		}
	}

	unreachable R;

	// pop from Qa and insert into Qb
	long a = 0, b = 1;
	while (!Q[a].empty()) {
		long v = Q[a].front();
		Q[a].pop();

		R.U.push_back(v);
		R.U_set.insert(v);
		for (long i = G.in_row[v]; i < G.in_row[v+1]; i++) {
			R.I_src.push_back(G.in_col[i]);
			R.I_dst.push_back(v);
		}

		long row_start = G.out_row[v], row_end = G.out_row[v+1];
		for (long i = G.out_row[v]; i < G.out_row[v+1]; i++) {
			long x = G.out_col[i];

			G.remove_edge(v, x);
			if (G.in_deg[x] == 0) {
				Q[b].push(x);
			}
		}

		// swap the queues
		if (Q[a].empty()) {
			long temp = a;
			a = b;
			b = a;
		}
	}

	return R;
}

// function to find the set of vertices which does not reach the sink w
unreachable find_unreachable_up(dyn_graph G, long ns, long *S, long w) {
	std::queue<long> Q[2];

	for (long i = 0; i < ns; i++) {
		long v = S[i];
		if (v != w) {
			// if v has no out-edges
			if (G.out_deg[v] == 0) {
				Q[0].push(v);
			}
		}
	}

	unreachable R;

	// pop from Qa and insert into Qb
	long a = 0, b = 1;
	while (!Q[a].empty()) {
		long v = Q[a].front();
		Q[a].pop();

		R.U.push_back(v);
		R.U_set.insert(v);
		for (long i = G.in_row[v]; i < G.in_row[v+1]; i++) {
			R.I_src.push_back(G.out_col[i]);
			R.I_dst.push_back(v);
		}

		for (long i = G.in_row[v]; i < G.in_row[v+1]; i++) {
			long x = G.in_col[i];

			G.remove_edge(x, v);
			if (G.out_deg[x] == 0) {
				Q[b].push(x);
			}
		}

		// swap the queues
		if (Q[a].empty()) {
			long temp = a;
			a = b;
			b = a;
		}
	}

	return R;
}
