# bellmanfordsp
Bellman Ford algorithm to find the shortest path in a digraph with negative-weight cycle detection
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

Enter the graph options
![BellmanFord_4](https://user-images.githubusercontent.com/117768679/223192308-f41c2054-49a7-4ff9-86eb-4907f387d569.PNG)

Enter the shortest path vertices:  Source Vertex and Target Vertex
![BellmanFord_3](https://user-images.githubusercontent.com/117768679/223192687-df44f513-9f9c-4e5d-84a5-cca6e308f41e.PNG)

Create a negative edge (weight = -1), shown in yellow
![BellmanFord_1](https://user-images.githubusercontent.com/117768679/223192882-b9bdf74a-b25c-4955-8354-9706a86e15fe.PNG)

Create a negative cycle by decreasing the weight (-18)
![BellmanFord_2](https://user-images.githubusercontent.com/117768679/223193350-acc98ef9-ed6a-420d-922f-45ab53d7321b.PNG)

Another Shortest Path
![BellmanFord_5](https://user-images.githubusercontent.com/117768679/223199517-ab282256-13ec-49e3-94b4-c26d91df1f0a.PNG)

Create a negative edge (weight = -2.00), shown in yellow
![BellmanFord_6](https://user-images.githubusercontent.com/117768679/223199782-cc040a71-c7cb-4eea-9f68-0c83e0a806fa.PNG)

Create a negative cycle by decreasing the weight (-3.00)
![BellmanFord_7](https://user-images.githubusercontent.com/117768679/223200157-181c5028-e9a2-4437-b938-9c224c0cd4bd.PNG)
