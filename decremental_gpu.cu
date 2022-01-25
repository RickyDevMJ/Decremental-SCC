#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/set_operations.h>
#include <thrust/execution_policy.h>

// todo : implement dyn_graph data-structure to represent graphs

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
	thrust::device_vector<long> U;
	thrust::device_vector<std::pair<long, long>> I;
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

// kernel corresponding to the helper function lift_up
__global__ void lift_up_kernel(scc_tree *T, unreachable &R) {

	long tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < T->num_children) {
		if (T->children[tid] == NULL)
				return;

		for (long j = 0; j < R.U.size(); j++) {
			if (T->children[tid] == t_array[R.U[j]]) {
				t_array[R.U[j]]->parent = T->parent;

				T->children[tid] = NULL;
				atomicSub(&T->num_children, 1);

				// todo : modify T.G

				if (T->parent != NULL)
					// todo : modify T->parent and T->parent.G

				break;
			}
		}
	}
}

// helper function to recursively lift up the deleted node in the SCC-tree
void lift_up(scc_tree *T, unreachable &R, scc_tree **t_array) {

	unreachable *d_R;
	cudaMalloc((void **)&d_R, sizeof(unreachable));
	cudaMemcpy(d_R, &R, sizeof(unreachable), cudaMemcpyHostToDevice);

	long nthreads = BLOCK_SIZE;
	long nblocks = (T->max_children - 1) / nthreads + 1;

	// todo: transfer tree T to device (d_T)
	lift_up_kernel<<< nthreads, nblocks >>>(d_T, d_R);

	// if node is not the root
	if (T->parent != NULL) {
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
	unreachable R1 = find_unreachable_down(T->G, 2, S, T->G.vertex_map[T->id]);
	unreachable R2 = find_unreachable_up(T->G, 2, S, T->G.vertex_map[T->id]);

	unreachable R;
	thrust::device_vector<long> U(R1.U.size() + R2.U.size());
	thrust::device_vector<std::pair<long, long>> I(R1.I.size() + R2.I.size());
	thrust::sort(thrust::device, R1.U.begin(), R1.U.end());
	thrust::sort(thrust::device, R2.U.begin(), R2.U.end());
	thrust::set_union(thrust::device, R1.U.begin(), R1.U.end(), R2.U.begin(), R2.U.end(), U.begin());
	thrust::sort(thrust::device, R1.I.begin(), R1.I.end());
	thrust::sort(thrust::device, R2.I.begin(), R2.I.end());
	thrust::set_union(thrust::device, R1.I.begin(), R1.I.end(), R2.I.begin(), R2.I.end(), I.begin());
	R.U = U;
	R.I = I;

	lift_up(T, R, t_array);
}

// function to find the set of vertices not reachable from the source w
unreachable find_unreachable_down(dyn_graph G, long ns, long *S, long w) {
	thrust::device_vector<long> Q[2];

	for (long i = 0; i < ns; i++) {
		long v = S[i];
		if (v != w) {
			// if v has no in-edges
			if (G.in_deg[v] == 0) {
				Q[0].push_back(v);
			}
		}
	}

	unreachable R;

	// pop from Qa and insert into Qb
	long a = 0, b = 1;
	while (!Q[a].empty()) {
		long v = Q[a][0];
		Q[a].erase(Q[a].begin());

		R.U.push_back(v);
		for (long i = G.out_row[v]; i < G.out_row[v+1]; i++) {
			R.I.push_back(v, G.out_col[i]));
		}

		long row_start = G.out_row[v], row_end = G.out_row[v+1];
		for (long i = G.out_row[v]; i < G.out_row[v+1]; i++) {
			long x = G.out_col[i];

			G.remove_edge(v, x);
			if (G.in_deg[x] == 0) {
				Q[b].push_back(x);
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
	thrust::device_vector<long> Q[2];

	for (long i = 0; i < ns; i++) {
		long v = S[i];
		if (v != w) {
			// if v has no out-edges
			if (G.out_deg[v] == 0) {
				Q[0].push_back(v);
			}
		}
	}

	unreachable R;

	// pop from Qa and insert into Qb
	long a = 0, b = 1;
	while (!Q[a].empty()) {
		long v = Q[a][0];
		Q[a].erase(Q[a].begin());

		R.U.push_back(v);
		for (long i = G.in_row[v]; i < G.in_row[v+1]; i++) {
			R.I.push_back(std::make_pair(G.in_col[i], v));
		}

		for (long i = G.in_row[v]; i < G.in_row[v+1]; i++) {
			long x = G.in_col[i];

			G.remove_edge(x, v);
			if (G.out_deg[x] == 0) {
				Q[b].push_back(x);
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
