#include <iostream>
#include "graph.h"

using namespace std;
typedef Graph<int,int,int> GraphType;

int main() {
  cout << "Running the README example" << endl;
	GraphType *g = new GraphType(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1);

	g -> add_node();
	g -> add_node();

	g -> add_tweights( 0,   /* capacities */  1, 5 );
	g -> add_tweights( 1,   /* capacities */  2, 6 );
	g -> add_edge( 0, 1,    /* capacities */  3, 4 );

	int flow = g -> maxflow();

  cout << "Flow " << flow << endl;
  cout << "Minimum cut:" << endl;
	if (g->what_segment(0) == GraphType::SOURCE)
		cout << "node0 is in the SOURCE set\n" << endl;
	else
		cout << "node0 is in the SINK set\n" << endl;
	if (g->what_segment(1) == GraphType::SOURCE)
		cout << "node1 is in the SOURCE set\n" << endl;
	else
		cout << "node1 is in the SINK set\n" << endl;

	delete g;
  g = 0;


  cout << "Running example on CLRS pp. 726 (3rd ed)" << endl;
  GraphType *c726 = new GraphType(4, 5);
  c726->add_node(5);
  c726->add_edge(1, 3, 12, 0);
  c726->add_edge(2, 1, 4, 0);
  c726->add_edge(3, 2, 9, 0);
  c726->add_edge(2, 4, 14, 0);
  c726->add_edge(4, 3, 7, 0);

  c726->add_tweights(1, 16, 0);
  c726->add_tweights(2, 13, 0);
  c726->add_tweights(3, 0, 20);
  c726->add_tweights(4, 0, 4);
  flow = c726->maxflow();
  cout << "Flow " << flow << endl;
  if (flow != 23) {
    cerr << "Error: expected flow 32" << endl;
  }
  delete c726;
  c726 = 0;

  cout << "Running example on CLRS pp. 728 (3rd ed)" << endl;
  GraphType *c728 = new GraphType(2, 1);
  c728->add_node(3);
  c728->add_edge(1, 2, 1, 0);
  const int mill = 1000000;
  c728->add_tweights(1, mill, mill);
  c728->add_tweights(1, mill, mill);
  flow = c728->maxflow();
  cout << "Flow " << flow << endl;
  if (flow != 2 * mill) {
    cerr << "Error: expected flow " << mill << endl;
  }
  delete c728;
  c728 = 0;

	return 0;
}


