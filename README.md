# bellmanfordsp
Bellman Ford algorithm to find the shortest path with negative-weight cycle detection
This is a web application written in Go and which uses the standard library package html/template to create the dynamic html for the web browser.
To use this program start the program with "go run sp.go", which starts the http server listening on localhost:8080.  Open a web browser and
address it to http://127.0.0.1:8080/graphoptions.  The user can select 2-500 vertices to create a randomly generated Euclidean graph in
the desired x-y boundary.  A minimum spanning tree is generated using the Prim algorithm and is plotted.  The user can then select a source
and target vertex for the shortest path (SP).  The Bellman Ford algorithm is used to determine the shortest path.  The vertices traversed and
the path's distance are displayed and plotted.  A negative cycle can be created along the shortest path by entering a vertex in Negative-edge
from and Negative-edge to and a negative weight.  The direction of the negative edge should be opposite the direction of the shortest path.
A negative cycle in a digraph is an ill-posed problem that has no solution.  It results in an infinite loop as the algorithm continuously finds
a shorter path from source to target.  The algorithm is terminated after the number of vertices in the path becomes the number of vertices in
the graph, as the shortest path cannot be any longer than this.  Note that the negative distance is highlighed in red and the Shortest Path 
Vertices are limited to the number of vertices in the graph.
