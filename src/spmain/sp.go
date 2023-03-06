/*
This is a web application.  The backend server is written in Go and uses the
html/package to create the html used by the web browser, which points to localhost:8080/bellmanfordsp.
Bellman Ford shortest paths finds the shortest paths (SP) between a source vertex and all other vertices.
Negative-weight edge cycles are detected.  The edges in the cycle are displayed.
Plot the SP showing the vertices and edges connecting the chosen source and target.
The user enters the following data in an html form:  #vertices and  x-y Euclidean bounds.
Euclidean graphs have no negative weights and no negative-weight cycles, so the user
can force negative-weight edge cycles by specifying the vertices and negative weight in the form.
A random number of vertices is chosen for the connection with a random start vertex.
The user can select the source and target vertices of the shortest path to find.  Their
coordinates are displayed as well as their distance.
*/

package main

import (
	"bufio"
	"container/heap"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"net/http"
	"os"
	"strconv"
	"strings"
	"text/template"
	"time"
)

const (
	addr                 = "127.0.0.1:8080"               // http server listen address
	fileBellmanFordSP    = "templates/bellmanfordsp.html" // html for bellmanford SP
	fileGraphOptions     = "templates/graphoptions.html"  // html for Graph Options
	patternBellmanFordSP = "/bellmanfordsp"               // http handler for bellmanford SP connections
	patternGraphOptions  = "/graphoptions"                // http handler for Graph Options
	rows                 = 300                            // #rows in grid
	columns              = rows                           // #columns in grid
	xlabels              = 11                             // # labels on x axis
	ylabels              = 11                             // # labels on y axis
	fileVerts            = "vertices.csv"                 // bounds and complex locations of vertices
)

// Edges are the vertices of the edge endpoints
type Edge struct {
	v int // one vertix
	w int // the other vertix
}

// Items are stored in the Priority Queue
type Item struct {
	Edge             // embedded field accessed with v,w
	index    int     // The index is used by Priority Queue update and is maintained by the heap.Interface
	distance float64 // Edge distance between vertices
}

// Priority Queue is a map of indexes and queue items and implements the heap.Interface
// A map is used instead of a slice so that it can be easily determined if an edge is in the queue
type PriorityQueue map[int]*Item

// Minimum spanning tree holds the edge vertices
type MST []*Edge

// Type to contain all the HTML template actions
type PlotT struct {
	Grid             []string // plotting grid
	Status           string   // status of the plot
	Xlabel           []string // x-axis labels
	Ylabel           []string // y-axis labels
	Distance         string   // Prim MST total distance (all the edges in MST)
	Vertices         string   // number of vertices
	Xmin             string   // x minimum endpoint in Euclidean graph
	Xmax             string   // x maximum endpoint in Euclidean graph
	Ymin             string   // y minimum endpoint in Euclidean graph
	Ymax             string   // y maximum endpoint in Euclidean graph
	StartLocation    string   // Prim MST start vertex location in x,y coordinates
	SourceLocation   string   // source vertex for bellmanford SP in x,y coordinates
	TargetLocation   string   // target or destination vertex for bellmanford SP in x,y coordinates
	Source           string   // source vertex for bellmanford SP 0-Vertices-1
	Target           string   // target vertex for bellmanford SP 0-Vertices-1
	DistanceSP       string   // shortest path distance (source->target)
	NegativeEdgeFrom string   // Negative edge from vertex
	NegativeEdgeTo   string   // Negative edge to vertex
	NegativeWeight   string   // Negative weight < 0
	PathSP           string   // List of vertices
	NegDistance      string   // Negative distance for SP
}

// Type to hold the minimum and maximum data values of the Euclidean graph
type Endpoints struct {
	xmin float64
	xmax float64
	ymin float64
	ymax float64
}

// PrimMST type for Minimum Spanning Tree methods
type PrimMST struct {
	graph      [][]float64  // matrix of vertices and their distance from each other
	location   []complex128 // complex point(x,y) coordinates of vertices
	mst        MST
	*Endpoints // Euclidean graph endpoints
	plot       *PlotT
}

// BellmanFordSP type for Shortest Path methods
type BellmanFordSP struct {
	edgeTo      []*Edge      // edge to vertex w
	distTo      []float64    // distance to w from source
	adj         [][]*Edge    // adjacency list
	mst         MST          // reference PrimMST
	graph       [][]float64  // reference PrimMST
	location    []complex128 // reference PrimMST
	plot        *PlotT       // reference PrimMST
	source      int          // start vertex for shortest path
	target      int          // end vertex for shortest path
	*Endpoints               // Euclidean graph endpoints
	negEdgeFrom int          // negative edge from vertex
	negEdgeTo   int          // negative edge to vertex
	negWeight   float64      // negative edge weight < 0
}

// global variables for parse and execution of the html template and MST construction
var (
	tmplForm *template.Template
)

// init parses the html template fileS
func init() {
	tmplForm = template.Must(template.ParseFiles(fileBellmanFordSP))
}

// generateVertices creates random vertices in the complex plane
func (p *PrimMST) generateVertices(r *http.Request) error {

	// if Source and Target have values, then graph was saved and
	// we are going to calculate the SP.
	sourceVert := r.PostFormValue("sourcevert")
	targetVert := r.PostFormValue("targetvert")
	if len(sourceVert) > 0 && len(targetVert) > 0 {
		f, err := os.Open(fileVerts)
		if err != nil {
			fmt.Printf("Open file %s error: %v\n", fileVerts, err)
		}
		defer f.Close()
		input := bufio.NewScanner(f)
		input.Scan()
		line := input.Text()
		// Each line has comma-separated values
		values := strings.Split(line, ",")
		var xmin, ymin, xmax, ymax float64
		if xmin, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			return err
		}

		if ymin, err = strconv.ParseFloat(values[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
			return err
		}
		if xmax, err = strconv.ParseFloat(values[2], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[2], err)
			return err
		}

		if ymax, err = strconv.ParseFloat(values[3], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[3], err)
			return err
		}
		p.Endpoints = &Endpoints{xmin: xmin, ymin: ymin, xmax: xmax, ymax: ymax}

		p.location = make([]complex128, 0)
		for input.Scan() {
			line := input.Text()
			// Each line has comma-separated values
			values := strings.Split(line, ",")
			var x, y float64
			if x, err = strconv.ParseFloat(values[0], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
				continue
			}
			if y, err = strconv.ParseFloat(values[1], 64); err != nil {
				fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
				continue
			}
			p.location = append(p.location, complex(x, y))
		}

		return nil
	}
	// Generate V vertices and locations randomly, get from HTML form
	// or read in from a previous graph when using a new start vertex.
	// Insert vertex complex coordinates into locations
	str := r.FormValue("xmin")
	xmin, err := strconv.ParseFloat(str, 64)
	if err != nil {
		fmt.Printf("String %s conversion to float error: %v\n", str, err)
		return err
	}

	str = r.FormValue("ymin")
	ymin, err := strconv.ParseFloat(str, 64)
	if err != nil {
		fmt.Printf("String %s conversion to float error: %v\n", str, err)
		return err
	}

	str = r.FormValue("xmax")
	xmax, err := strconv.ParseFloat(str, 64)
	if err != nil {
		fmt.Printf("String %s conversion to float error: %v\n", str, err)
		return err
	}

	str = r.FormValue("ymax")
	ymax, err := strconv.ParseFloat(str, 64)
	if err != nil {
		fmt.Printf("String %s conversion to float error: %v\n", str, err)
		return err
	}

	// Check if xmin < xmax and ymin < ymax and correct if necessary
	if xmin >= xmax {
		xmin, xmax = xmax, xmin
	}
	if ymin >= ymax {
		ymin, ymax = ymax, ymin
	}

	p.Endpoints = &Endpoints{xmin: xmin, ymin: ymin, xmax: xmax, ymax: ymax}

	vertices := r.FormValue("vertices")
	verts, err := strconv.Atoi(vertices)
	if err != nil {
		fmt.Printf("String %s conversion to int error: %v\n", vertices, err)
		return err
	}

	delx := xmax - xmin
	dely := ymax - ymin
	// Generate vertices
	p.location = make([]complex128, verts)
	for i := 0; i < verts; i++ {
		x := xmin + delx*rand.Float64()
		y := ymin + dely*rand.Float64()
		p.location[i] = complex(x, y)
	}

	// Save the endpoints and vertex locations to a csv file
	f, err := os.Create(fileVerts)
	if err != nil {
		fmt.Printf("Create file %s error: %v\n", fileVerts, err)
		return err
	}
	defer f.Close()
	// Save the endpoints
	fmt.Fprintf(f, "%f,%f,%f,%f\n", p.xmin, p.ymin, p.xmax, p.ymax)
	// Save the vertex locations as x,y
	for _, z := range p.location {
		fmt.Fprintf(f, "%f,%f\n", real(z), imag(z))
	}

	return nil
}

// findDistances find distances between vertices and insert into graph
func (p *PrimMST) findDistances() error {

	verts := len(p.location)
	// Store distances between vertices for Euclidean graph
	p.graph = make([][]float64, verts)
	for i := 0; i < verts; i++ {
		p.graph[i] = make([]float64, verts)
	}

	for i := 0; i < verts; i++ {
		for j := i + 1; j < verts; j++ {
			distance := cmplx.Abs(p.location[i] - p.location[j])
			p.graph[i][j] = distance
			p.graph[j][i] = distance
		}
	}
	for i := 0; i < verts; i++ {
		p.graph[i][i] = math.MaxFloat64
	}

	return nil
}

// A PriorityQueue implements heap.Interface and holds Items
// Len returns length of queue.
func (pq PriorityQueue) Len() int {
	return len(pq)
}

// Less returns Item weight[i] less than Item weight[j]
func (pq PriorityQueue) Less(i, j int) bool {
	return pq[i].distance < pq[j].distance
}

// Swap swaps Item[i] and Item[j]
func (pq PriorityQueue) Swap(i, j int) {
	pq[i], (pq)[j] = pq[j], pq[i]
	pq[i].index = i
	pq[j].index = j
}

// Push inserts an Item in the queue
func (pq *PriorityQueue) Push(x interface{}) {
	n := len(*pq)
	item := x.(*Item)
	item.index = n
	(*pq)[n] = item
}

// Pop removes an Item from the queue and returns it
func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	old[n-1] = nil
	item.index = -1
	delete(*pq, n-1)
	return item
}

// update modifies the distance and value of an Item in the queue
func (pq *PriorityQueue) update(item *Item, distance float64) {
	item.distance = distance
	heap.Fix(pq, item.index)
}

// findMST finds the minimum spanning tree (MST) using Prim's algorithm
func (p *PrimMST) findMST() error {
	vertices := len(p.location)
	p.mst = make(MST, vertices)
	marked := make([]bool, vertices)
	distTo := make([]float64, vertices)
	for i := range distTo {
		distTo[i] = math.MaxFloat64
	}
	// Create a priority queue, put the items in it, and establish
	// the priority queue (heap) invariants.
	pq := make(PriorityQueue)

	visit := func(v int) {
		marked[v] = true
		// find shortest distance from vertex v to w
		for w, dist := range p.graph[v] {
			// Check if already in the MST
			if marked[w] {
				continue
			}
			if dist < distTo[w] {
				// Edge to w is new best connection from MST to w
				p.mst[w] = &Edge{v: v, w: w}
				distTo[w] = dist
				// Check if already in the queue and update
				item, ok := pq[w]
				// update
				if ok {
					pq.update(item, dist)
				} else { // insert
					item = &Item{Edge: Edge{v: v, w: w}, distance: dist}
					heap.Push(&pq, item)
				}
			}
		}
	}

	// Starting index is 0, distance is MaxFloat64, put it in the queue
	distTo[0] = math.MaxFloat64
	pq[0] = &Item{index: 0, distance: math.MaxFloat64, Edge: Edge{v: 0, w: 0}}
	heap.Init(&pq)

	// Loop until the queue is empty and the MST is finished
	for pq.Len() > 0 {
		item := heap.Pop(&pq).(*Item)
		visit(item.w)
	}

	return nil
}

// plotMST draws the MST onto the grid
func (p *PrimMST) plotMST(status []string) error {

	// Apply the parsed HTML template to plot object
	// Construct x-axis labels, y-axis labels, status message

	var (
		xscale   float64
		yscale   float64
		distance float64
	)
	p.plot = &PlotT{}
	p.plot.Grid = make([]string, rows*columns)
	p.plot.Xlabel = make([]string, xlabels)
	p.plot.Ylabel = make([]string, ylabels)

	// Calculate scale factors for x and y
	xscale = (columns - 1) / (p.xmax - p.xmin)
	yscale = (rows - 1) / (p.ymax - p.ymin)

	// Insert the mst vertices and edges in the grid
	// loop over the MST vertices

	// color the vertices black
	// color the edges connecting the vertices gray
	// color the MST start vertex green
	// create the line y = mx + b for each edge
	// translate complex coordinates to row/col on the grid
	// translate row/col to slice data object []string Grid
	// CSS selectors for background-color are "vertex", "startvertexMSS", and "edge"

	beginEP := complex(p.xmin, p.ymin)  // beginning of the Euclidean graph
	endEP := complex(p.xmax, p.ymax)    // end of the Euclidean graph
	lenEP := cmplx.Abs(endEP - beginEP) // length of the Euclidean graph

	for _, e := range p.mst[1:] {

		// Insert the edge between the vertices v, w.  Do this before marking the vertices.
		// CSS colors the edge gray.
		beginEdge := p.location[e.v]
		endEdge := p.location[e.w]
		lenEdge := cmplx.Abs(endEdge - beginEdge)
		distance += lenEdge
		ncells := int(columns * lenEdge / lenEP) // number of points to plot in the edge

		beginX := real(beginEdge)
		endX := real(endEdge)
		deltaX := endX - beginX
		stepX := deltaX / float64(ncells)

		beginY := imag(beginEdge)
		endY := imag(endEdge)
		deltaY := endY - beginY
		stepY := deltaY / float64(ncells)

		// loop to draw the edge
		x := beginX
		y := beginY
		for i := 0; i < ncells; i++ {
			row := int((p.ymax-y)*yscale + .5)
			col := int((x-p.xmin)*xscale + .5)
			p.plot.Grid[row*columns+col] = "edge"
			x += stepX
			y += stepY
		}

		// Mark the edge start vertex v.  CSS colors the vertex black.
		row := int((p.ymax-beginY)*yscale + .5)
		col := int((beginX-p.xmin)*xscale + .5)
		p.plot.Grid[row*columns+col] = "vertex"

		// Mark the edge end vertex w.  CSS colors the vertex black.
		row = int((p.ymax-endY)*yscale + .5)
		col = int((endX-p.xmin)*xscale + .5)
		p.plot.Grid[row*columns+col] = "vertex"
	}

	// Mark the MST start vertex.  CSS colors the vertex green.
	x := real(p.location[0])
	y := imag(p.location[0])
	p.plot.StartLocation = fmt.Sprintf("(%.2f, %.2f)", x, y)
	row := int((p.ymax-y)*yscale + .5)
	col := int((x-p.xmin)*xscale + .5)
	p.plot.Grid[row*columns+col] = "startvertexMSS"
	p.plot.Grid[(row+1)*columns+col] = "startvertexMSS"
	p.plot.Grid[(row-1)*columns+col] = "startvertexMSS"
	p.plot.Grid[row*columns+col+1] = "startvertexMSS"
	p.plot.Grid[row*columns+col-1] = "startvertexMSS"

	// Construct x-axis labels
	incr := (p.xmax - p.xmin) / (xlabels - 1)
	x = p.xmin
	// First label is empty for alignment purposes
	for i := range p.plot.Xlabel {
		p.plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (p.ymax - p.ymin) / (ylabels - 1)
	y = p.ymin
	for i := range p.plot.Ylabel {
		p.plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Distance of the MST
	p.plot.Distance = fmt.Sprintf("%.2f", distance)

	// Endpoints and Vertices
	p.plot.Vertices = strconv.Itoa(len(p.location))
	p.plot.Xmin = fmt.Sprintf("%.2f", p.xmin)
	p.plot.Xmax = fmt.Sprintf("%.2f", p.xmax)
	p.plot.Ymin = fmt.Sprintf("%.2f", p.ymin)
	p.plot.Ymax = fmt.Sprintf("%.2f", p.ymax)

	return nil
}

// findSP constructs the shortest path from source to target
func (bfsp *BellmanFordSP) findSP(r *http.Request) error {
	// need both source and target vertices for the shortest path
	sourceVert := r.PostFormValue("sourcevert")
	targetVert := r.PostFormValue("targetvert")
	var err error
	if len(sourceVert) == 0 || len(targetVert) == 0 {
		return fmt.Errorf("source and/or target vertices not set")
	}
	bfsp.source, err = strconv.Atoi(sourceVert)
	if err != nil {
		fmt.Printf("source vertex Atoi error: %v\n", err)
		return err
	}
	bfsp.target, err = strconv.Atoi(targetVert)
	if err != nil {
		fmt.Printf("target vertex Atoi error: %v\n", err)
		return err
	}

	vertices := len(bfsp.location)
	if bfsp.source == bfsp.target || bfsp.source < 0 || bfsp.target < 0 ||
		bfsp.source > vertices-1 || bfsp.target > vertices-1 {
		return fmt.Errorf("source and/or target vertices are invalid")
	}

	bfsp.edgeTo = make([]*Edge, vertices)
	bfsp.distTo = make([]float64, vertices)
	for i := range bfsp.distTo {
		bfsp.distTo[i] = math.MaxFloat64
	}
	// Starting index is source, distance to itself is 0
	bfsp.distTo[bfsp.source] = 0.0

	// Create the adjacency list
	bfsp.adj = make([][]*Edge, vertices)
	for i := range bfsp.adj {
		bfsp.adj[i] = make([]*Edge, 0)
	}
	for _, e := range bfsp.mst[1:] {
		bfsp.adj[e.v] = append(bfsp.adj[e.v], e)
		bfsp.adj[e.w] = append(bfsp.adj[e.w], e)
	}

	// relax the vertex and obtain the shortest path if no negative-weight cycles
	relax := func(v int) {
		// find shortest distance from source to w
		for _, e := range bfsp.adj[v] {
			// Determine v and w on the edge
			w := e.w
			if e.w == v {
				w = e.v
				e.v, e.w = e.w, e.v
			}

			newDistance := bfsp.distTo[v] + bfsp.graph[v][w]
			if bfsp.distTo[w] > newDistance {
				// Edge to w is new best connection from source to w
				bfsp.edgeTo[w] = e
				bfsp.distTo[w] = newDistance
			}
		}
	}

	// detect negative-weight cycle and return edge
	detect := func(v int) []*Edge {
		edges := make([]*Edge, 0)
		for _, e := range bfsp.adj[v] {
			// Determine v and w on the edge
			w := e.w
			if e.w == v {
				w = e.v
				e.v, e.w = e.w, e.v
			}

			// if distance decreases, this indicates a negative-weight cycle
			newDistance := bfsp.distTo[v] + bfsp.graph[v][w]
			if bfsp.distTo[w] > newDistance {
				edges = append(edges, e)
			}
		}
		return edges
	}

	// Insert the negative-weight if the HTML form controls are set
	negEdgeFrom := r.PostFormValue("negedgefrom")
	negEdgeTo := r.PostFormValue("negedgeto")
	negWeight := r.PostFormValue("negweight")
	if len(negEdgeFrom) > 0 && len(negEdgeTo) > 0 && len(negWeight) > 0 &&
		negEdgeFrom != "0" && negEdgeTo != "0" && negWeight != "0.00" {
		bfsp.negEdgeFrom, err = strconv.Atoi(negEdgeFrom)
		if err != nil {
			fmt.Printf("negative-edge from Atoi error: %v\n", err)
			return err
		}
		bfsp.negEdgeTo, err = strconv.Atoi(negEdgeTo)
		if err != nil {
			fmt.Printf("negative-edge to Atoi error: %v\n", err)
			return err
		}
		bfsp.negWeight, err = strconv.ParseFloat(negWeight, 64)
		if err != nil {
			fmt.Printf("negative weight ParseFloat error: %v\n", err)
			return err
		}

		if bfsp.negEdgeFrom == bfsp.negEdgeTo || bfsp.negEdgeFrom < 0 || bfsp.negEdgeTo < 0 ||
			bfsp.negEdgeFrom > vertices-1 || bfsp.negEdgeTo > vertices-1 {
			return fmt.Errorf("negative-edge from/to vertices are invalid")
		}
		if bfsp.negWeight >= 0 {
			return fmt.Errorf("negative-edge weight must be less than zero")
		}

		edge := &Edge{v: bfsp.negEdgeFrom, w: bfsp.negEdgeTo}
		bfsp.adj[bfsp.negEdgeFrom] = append(bfsp.adj[bfsp.negEdgeFrom], edge)
		bfsp.graph[bfsp.negEdgeFrom][bfsp.negEdgeTo] = bfsp.negWeight

	}

	// Perform the number of vertices minus one passes over all the vertices
	// Finds the shortest path if no negative-weight cycles are present
	for pass := 0; pass < vertices-1; pass++ {
		for v := 0; v < vertices; v++ {
			relax(v)
		}
	}

	// Detect negative-weight cycle that is reachable from the source vertex
	cycles := make([]*Edge, 0)
	for v := 0; v < vertices; v++ {
		edges := detect(v)
		if len(edges) > 0 {
			cycles = append(cycles, edges...)
		}
	}
	// convert edges to error
	errors := make([]string, 0)
	for _, e := range cycles {
		errors = append(errors, fmt.Sprintf("(%d,%d)", e.v, e.w))
	}
	if len(errors) > 0 {
		return fmt.Errorf("negative-weight cycles in edges: " + strings.Join(errors, ", "))
	}

	return nil
}

// plotSP draws the shortest path from source to target in the grid and any
// edges in negative-weight cycles
func (bfsp *BellmanFordSP) plotSP() error {
	// check if the target was found in findSP
	if len(bfsp.distTo) == 0 || bfsp.distTo[bfsp.target] == math.MaxFloat64 {
		return fmt.Errorf("distance to vertex %d not found", bfsp.target)
	}

	var (
		distance float64 = 0.0
		edges    []*Edge = make([]*Edge, 0)
	)

	// Calculate scale factors for x and y
	xscale := (columns - 1) / (bfsp.xmax - bfsp.xmin)
	yscale := (rows - 1) / (bfsp.ymax - bfsp.ymin)

	beginEP := complex(bfsp.xmin, bfsp.ymin) // beginning of the Euclidean graph
	endEP := complex(bfsp.xmax, bfsp.ymax)   // end of the Euclidean graph
	lenEP := cmplx.Abs(endEP - beginEP)      // length of the Euclidean graph

	e := bfsp.edgeTo[bfsp.target]
	if e.w != bfsp.target {
		e.v, e.w = e.w, e.v
	}
	// start at the target and loop until source vertex is plotted to the grid
	nedges := 0
	for {
		v := e.v
		w := e.w

		edges = append(edges, e)

		start := bfsp.location[v]
		end := bfsp.location[w]
		x1 := real(start)
		y1 := imag(start)
		x2 := real(end)
		y2 := imag(end)
		lenEdge := cmplx.Abs(end - start)
		distance += bfsp.graph[v][w]
		ncells := int(columns * lenEdge / lenEP) // number of points to plot in the edge

		deltaX := x2 - x1
		stepX := deltaX / float64(ncells)

		deltaY := y2 - y1
		stepY := deltaY / float64(ncells)

		// loop to draw the SP edge; CSS colors the edge Orange
		x := x1
		y := y1
		for i := 0; i < ncells; i++ {
			row := int((bfsp.ymax-y)*yscale + .5)
			col := int((x-bfsp.xmin)*xscale + .5)
			bfsp.plot.Grid[row*columns+col] = "edgeSP"
			x += stepX
			y += stepY
		}

		// Mark the edge start vertex v.  CSS colors the vertex Black.
		row := int((bfsp.ymax-y1)*yscale + .5)
		col := int((x1-bfsp.xmin)*xscale + .5)
		bfsp.plot.Grid[row*columns+col] = "vertex"

		// Mark the edge end vertex w.  CSS colors the vertex Black.
		row = int((bfsp.ymax-y2)*yscale + .5)
		col = int((x2-bfsp.xmin)*xscale + .5)
		bfsp.plot.Grid[row*columns+col] = "vertex"

		vertices := len(bfsp.location)
		nedges++

		// Exit the loop if source is reached, we have the SP.
		// Or an infinite loop, stop after the number of vertices
		// in the graph is acquired.
		if e.v == bfsp.source || nedges == vertices-1 {
			break
		}

		// move forward to the next edge
		e = bfsp.edgeTo[v]
		if e.w != v {
			e.v, e.w = e.w, e.v
		}
	}

	// Draw the negative-weight edge
	v := bfsp.negEdgeFrom
	w := bfsp.negEdgeTo
	start := bfsp.location[v]
	end := bfsp.location[w]
	x1 := real(start)
	y1 := imag(start)
	x2 := real(end)
	y2 := imag(end)
	lenEdge := cmplx.Abs(end - start)
	ncells := int(columns * lenEdge / lenEP) // number of points to plot in the edge

	deltaX := x2 - x1
	stepX := deltaX / float64(ncells)

	deltaY := y2 - y1
	stepY := deltaY / float64(ncells)

	// loop to draw the edge; CSS colors the cycle edge Yellow
	x := x1
	y := y1
	for i := 0; i < ncells; i++ {
		row := int((bfsp.ymax-y)*yscale + .5)
		col := int((x-bfsp.xmin)*xscale + .5)
		bfsp.plot.Grid[row*columns+col] = "edgeCycle"
		x += stepX
		y += stepY
	}

	// Mark the end vertices of the shortest path
	e = bfsp.edgeTo[bfsp.target]
	x = real(bfsp.location[e.w])
	y = imag(bfsp.location[e.w])
	// Mark the SP end vertex.  CSS colors the vertex Red.
	row := int((bfsp.ymax-y)*yscale + .5)
	col := int((x-bfsp.xmin)*xscale + .5)
	bfsp.plot.Grid[row*columns+col] = "vertexSP2"
	bfsp.plot.Grid[(row+1)*columns+col] = "vertexSP2"
	bfsp.plot.Grid[(row-1)*columns+col] = "vertexSP2"
	bfsp.plot.Grid[row*columns+col+1] = "vertexSP2"
	bfsp.plot.Grid[row*columns+col-1] = "vertexSP2"

	bfsp.plot.TargetLocation = fmt.Sprintf("(%.2f, %.2f)", x, y)
	bfsp.plot.Target = strconv.Itoa(e.w)

	// Mark the SP start vertex.  CSS colors the vertex Blue.
	x = real(bfsp.location[bfsp.source])
	y = imag(bfsp.location[bfsp.source])
	row = int((bfsp.ymax-y)*yscale + .5)
	col = int((x-bfsp.xmin)*xscale + .5)
	bfsp.plot.Grid[row*columns+col] = "vertexSP1"
	bfsp.plot.Grid[(row+1)*columns+col] = "vertexSP1"
	bfsp.plot.Grid[(row-1)*columns+col] = "vertexSP1"
	bfsp.plot.Grid[row*columns+col+1] = "vertexSP1"
	bfsp.plot.Grid[row*columns+col-1] = "vertexSP1"

	bfsp.plot.SourceLocation = fmt.Sprintf("(%.2f, %.2f)", x, y)
	bfsp.plot.Source = strconv.Itoa(bfsp.source)

	// Distance of the SP
	bfsp.plot.DistanceSP = fmt.Sprintf("%.2f", distance)
	if distance < 0.0 {
		bfsp.plot.NegDistance = "negativedistance"
	}

	// Negative-edge vertices and weight
	bfsp.plot.NegativeEdgeFrom = strconv.Itoa(bfsp.negEdgeFrom)
	bfsp.plot.NegativeEdgeTo = strconv.Itoa(bfsp.negEdgeTo)
	bfsp.plot.NegativeWeight = fmt.Sprintf("%.2f", bfsp.negWeight)

	// list of vertices on the Shortest Path
	n := len(edges)
	verts := make([]string, n+1)
	for i := 0; i < n; i++ {
		verts[n-i] = strconv.Itoa(edges[i].w)
	}
	verts[0] = strconv.Itoa(edges[n-1].v)

	bfsp.plot.PathSP = strings.Join(verts, ", ")

	return nil

}

// HTTP handler for /graphoptions connections
func handleGraphOptions(w http.ResponseWriter, r *http.Request) {
	http.ServeFile(w, r, "templates/graphoptions.html")
}

// HTTP handler for /bellmanford connections
func handleBellmanFordSP(w http.ResponseWriter, r *http.Request) {

	// Create the Prim MST instance
	primmst := &PrimMST{}

	// Create the Bellman Ford SP instance
	bellmanfordsp := &BellmanFordSP{}

	// Accumulate error
	status := make([]string, 0)

	// Generate V vertices and locations randomly, get from HTML form
	// or read in from a previous graph when using a new start vertex.
	// Insert vertex complex coordinates into locations
	err := primmst.generateVertices(r)
	if err != nil {
		fmt.Printf("generateVertices error: %v\n", err)
		status = append(status, err.Error())
	}

	// Insert distances into graph
	err = primmst.findDistances()
	if err != nil {
		fmt.Printf("findDistances error: %v", err)
		status = append(status, err.Error())
	}

	// Find MST and save in PrimMST.mst
	err = primmst.findMST()
	if err != nil {
		fmt.Printf("findMST error: %v\n", err)
		status = append(status, err.Error())
	}

	// Assign vertex locations to bellmanfordsp so it can use x,y coordinates of vertices
	bellmanfordsp.location = primmst.location
	// Assign graph to bellmanfordsp so it can use distances between vertices
	bellmanfordsp.graph = primmst.graph
	// Assign MST to bellmanfordsp so it can use it to construct adj
	bellmanfordsp.mst = primmst.mst
	// Assign endpoints to bellmanfordsp for plotting on the grid
	bellmanfordsp.Endpoints = primmst.Endpoints

	// Find the Shortest Path
	err = bellmanfordsp.findSP(r)
	if err != nil {
		fmt.Printf("findSP error: %v\n", err)
		status = append(status, err.Error())
	}

	// Draw MST into 300 x 300 cell 2px grid
	// Construct x-axis labels, y-axis labels, status message
	err = primmst.plotMST(status)
	if err != nil {
		fmt.Printf("plotMST error: %v\n", err)
		status = append(status, err.Error())
	}

	// Assign plot to bellmanfordsp
	bellmanfordsp.plot = primmst.plot

	// Draw SP into 300 x 300 cell 2px grid
	err = bellmanfordsp.plotSP()
	if err != nil {
		fmt.Printf("plotSP error: %v\n", err)
		status = append(status, err.Error())
	}

	// Status
	if len(status) > 0 {
		bellmanfordsp.plot.Status = strings.Join(status, ", ")
	} else {
		bellmanfordsp.plot.Status = "Enter Source and Target Vertices (0-V-1) for another SP"
	}

	// Write to HTTP using template and grid
	if err := tmplForm.Execute(w, primmst.plot); err != nil {
		log.Fatalf("Write to HTTP output using template with grid error: %v\n", err)
	}
}

// main sets up the http handlers, listens, and serves http clients
func main() {
	rand.Seed(time.Now().Unix())
	// Set up http servers with handler for Graph Options and Bellman Ford SP
	http.HandleFunc(patternBellmanFordSP, handleBellmanFordSP)
	http.HandleFunc(patternGraphOptions, handleGraphOptions)
	fmt.Printf("Bellman Ford Shortest Path Server listening on %v.\n", addr)
	http.ListenAndServe(addr, nil)
}
