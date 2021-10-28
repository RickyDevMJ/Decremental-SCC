#include <vector>
#include <queue>

class dyn_graph {

	public:

		long n;
		long *in_row; 
		long *out_row;
		long *in_col; 
		long *out_col;
		long *in_deg;
		long *out_deg;

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
};

void find_unreachable_down(dyn_graph G, long ns, long *S, long w) {
    std::queue<long> Q;

    for (long i = 0; i < ns; i++) {
		long v = S[i];
        if (v != w) {
			// if v has no in-edges
			if (G.in_deg[v] == 0) {
				Q.push(v);
			}
		}
    }

	// vector is used instead of set since input graph is DAG
	std::vector<long> U, I_src, I_dst;

	while (!Q.empty()) {
		long v = Q.front();
		Q.pop();

		U.push_back(v);
		for (long i = G.in_row[v]; i < G.in_row[v+1]; i++) {
			I_src.push_back(G.in_col[i]);
			I_dst.push_back(v);
		}

		long row_start = G.out_row[v], row_end = G.out_row[v+1];
		for (long i = G.out_row[v]; i < G.out_row[v+1]; i++) {
			long x = G.out_col[i];

			G.remove_edge(v, x);
			if (G.in_deg[x] == 0) {
				Q.push(x);
			}
		}
	}

	// return U,I
}