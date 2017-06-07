/*
This is a Matlab interface to the maxflow/mincut algorithm defined in graph.h:

	This software library implements the maxflow algorithm
	described in

		An Experimental Comparison of Min-Cut/Max-Flow Algorithms
		for Energy Minimization in Vision.
		Yuri Boykov and Vladimir Kolmogorov.
		In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 
		September 2004

	This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
	at Siemens Corporate Research. To make it available for public use,
	it was later reimplemented by Vladimir Kolmogorov based on open publications.

	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.

USAGE:
See .m file.

BUGS:
memory management. Should somehow use Matlab's memory allocation instead of default new/delete.
	
*/  

#include "mex.h"
#include "graph.h"
#include <stdio.h>

/* 

the function will receive 2 arguments: edges tweights
edges is a matrix that specifies the edges. each row [e1 e2 c1 c2]
establishes an edge from e1 to e2 with capacity c1 
and from e2 to e1 with capacity c2
tweights ("terminal weights" I guess) specifies the edges from/to the
terminal nodes. the row [e s t] sets the capacity from S to e and from e to t
to s and t respectively. each edge should appear at most one time in this
matrix.
 
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //printf("DEBUG: Hello and welcome.\n");
	
	if (nrhs != 2) {
		mexErrMsgTxt("Expecting two arguments.");	
	}
	const mxArray *edges = prhs[0];
	// begin by verifying that the arguments are ok, since mexErrMsgTxt()
	// does not free the memory allocated with new
	// (actually this is a problem anyway since the function is stopped
	// if the user presses Ctrl-C, so we'll have to deal with it)
 	if (!mxIsDouble(edges)) {
		mexErrMsgTxt("EDGES should be of class double.");
	}
	int edges_m = mxGetM(edges);
	if (mxGetN(edges) != 4 || edges_m == 0) {
		mexErrMsgTxt("EDGES should be a N x 4 non-empty matrix.");
	} 
	int max_v = 0; // the maximum node number seen thus far
	double *edges_p = mxGetPr(edges);
	int i;
	for (i=0; i<edges_m; i++) {
		int e1 = (int)edges_p[i];
		int e2 = (int)edges_p[i+edges_m];
		if (e1<1 || e2<1) {
			mexErrMsgTxt("Vertex numbers must be positive integers.");
		}
		if (e1 > max_v) max_v = e1;
		if (e2 > max_v) max_v = e2;
	}	

	// check the TWEIGHTS matrix
	const mxArray *tweights = prhs[1];
	double *tweights_p = mxGetPr(tweights);
	if (!mxIsDouble(tweights)) {
		mexErrMsgTxt("TWEIGHTS should be of class double.");
	}
        int tweights_m = mxGetM(tweights);
        if (mxGetN(tweights) != 3 || tweights_m == 0) {
                mexErrMsgTxt("TWEIGHTS should be a N x 3 non-empty matrix.");
        }
	for (i=0; i<tweights_m; i++) {
		int v = (int)tweights_p[i];
		if (v<1) {
			mexErrMsgTxt("Vertex numbers must be positive integers.");
		}
		if (v > max_v) max_v = v;
	}

	// ok, done checking, now do it
	//printf("Number of vertices: %d\n", max_v);
	//printf("Number of edges: %d\n", edges_m);
	//printf("Number of tweights: %d\n", tweights_m);
	
	Graph *graph = new Graph();
	// mapping from numbers to vertex ids
	// node_array[i] is the node_id of node number i+1
	Graph::node_id *node_array = new Graph::node_id[max_v];
	
	for (i=0; i<max_v; i++) {
		node_array[i] = graph->add_node();
	}

	// add edges
	for (i=0; i<edges_m; i++) {
		int v1 = (int)edges_p[i];
		int v2 = (int)edges_p[i+edges_m];
		double w1 = edges_p[i+2*edges_m];
		double w2 = edges_p[i+3*edges_m];
		graph->add_edge(node_array[v1-1], node_array[v2-1], w1, w2);
		//printf("Added edges (%d, %d, %g, %g)\n", v1, v2, w1, w2);
		//fflush(stdout);
	}

	//printf("Done adding edges.\n");
	
	// add edges to terminal nodes
	for (i=0; i<tweights_m; i++) {
		// ASSUME no vertex occurs twice.
		int v = (int)tweights_p[i];
		double s = tweights_p[i+tweights_m];
		double t = tweights_p[i+2*tweights_m];
		graph->set_tweights(node_array[v-1], s, t);
		//printf("Set tweights (%d, %g, %g)\n", v, s, t);
	}


	// really do it
	Graph::flowtype maxflow = graph->maxflow();
	
	// put result in output parameters (if there are any)
	// put maxflow in 1st parameter
	if (nlhs > 0) {
		plhs[0] = mxCreateDoubleScalar(maxflow);
	}
	if (nlhs > 1) {
		// put labeling in 2nd parameter
		plhs[1] = mxCreateDoubleMatrix(1, max_v, mxREAL);
		double *p = mxGetPr(plhs[1]);
		for (i=0; i<max_v; i++) {
			p[i] = (graph->what_segment(node_array[i]) == Graph::SOURCE)?1:2;
		}
	}
	
	// clean up
	delete graph;
	delete [] node_array;	
}

