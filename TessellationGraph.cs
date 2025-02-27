using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

// Rhino namespaces for geometry and UI support.
using Rhino;
using Rhino.Geometry;
using Rhino.UI;

namespace Liz
{
    /// <summary>
    /// TessellationGraph builds a planar graph with a tessellated structure.
    /// Each node in the graph has a fixed number of connections (degree),
    /// and each face (polygon) in the tessellation has a fixed number of edges.
    /// </summary>
    public class TessellationGraph
    {
        // -------------------- Private Members --------------------

        /// <summary>
        /// List of all nodes (vertices) in the graph.
        /// </summary>
        private List<Node> _nodes;

        /// <summary>
        /// List of polygons (faces) in the graph, where each polygon is a list of node indices.
        /// </summary>
        private List<List<int>> _polygons;

        /// <summary>
        /// The number of sides for the initial (inner) polygon.
        /// </summary>
        private int _innerPolygon;

        /// <summary>
        /// The number of edges each polygon face should have.
        /// </summary>
        private int _outerPolygon;

        /// <summary>
        /// The default degree (number of connections) for each node.
        /// </summary>
        private int _defaultGraphDegree;

        /// <summary>
        /// The radius used for positioning nodes in a circular layout.
        /// </summary>
        private double _radius = 100.0;

        /// <summary>
        /// The current number of levels (layers) in the graph.
        /// </summary>
        private int _endLevel = 1;

        /// <summary>
        /// The starting index of nodes that have not been updated (newly added).
        /// </summary>
        private int _pendingNodesStart;

        /// <summary>
        /// The ending index (exclusive) of nodes that have not been updated.
        /// </summary>
        private int _pendingNodesEnd;

        /// <summary>
        /// Optional list to specify special graph degrees for specific levels.
        /// If a level is not specified, DefaultGraphDegree is used.
        /// </summary>
        private List<int> _specialGraphDegrees;

        // -------------------- Public Properties --------------------

        /// <summary>
        /// Returns the maximum level (layer) index in the graph.
        /// </summary>
        public int MaxLevel
        {
            get { return _endLevel - 1; }
        }

        // -------------------- Constructors --------------------

        /// <summary>
        /// Main constructor for TessellationGraph.
        /// Initializes the graph using an initial inner polygon and optional special degree settings.
        /// </summary>
        /// <param name="innerPolygon">Number of nodes/sides in the initial polygon.</param>
        /// <param name="outerPolygon">Number of edges for each tessellated face.</param>
        /// <param name="defaultGraphDegree">Default degree (number of connections) for each node.</param>
        /// <param name="specialGraphDegrees">
        /// Optional list of tuples indicating (level, degree) for levels that need special treatment.
        /// </param>
        public TessellationGraph(int innerPolygon, int outerPolygon, int defaultGraphDegree, List<Tuple<int, int>> specialGraphDegrees = null)
        {
            // Set basic parameters.
            _innerPolygon = innerPolygon;
            _outerPolygon = outerPolygon;
            _defaultGraphDegree = defaultGraphDegree; // Later: can validate that degree is >= 4.

            // Initialize node and polygon lists.
            _nodes = new List<Node>();
            _polygons = new List<List<int>>();
            _polygons.Add(new List<int>()); // Start with one polygon (the inner polygon).

            // Create the inner polygon by instantiating nodes and setting up circular neighbor indices.
            for (int i = 0; i < _innerPolygon; i++)
            {
                // Determine the index for the previous (before) and next node in the loop.
                int before = (i != 0) ? i - 1 : _innerPolygon - 1;
                int next = (i + 1) % _innerPolygon;
                // Create a new node for the first layer.
                var node = new Node(before, next);
                // Add the node index to the inner polygon.
                _polygons[0].Add(i);
                // Add the node to the node list.
                _nodes.Add(node);
            }

            // Connect each node to its next neighbor to form the closed loop of the inner polygon.
            for (int i = 0; i < _innerPolygon; i++)
            {
                Connect(i, (i + 1) % _innerPolygon);
            }

            // Handle special graph degree settings, if provided.
            if (specialGraphDegrees != null)
            {
                int maxDegreeIndex = 0;
                // Find the maximum level index specified.
                for (int i = 0; i < specialGraphDegrees.Count; i++)
                {
                    maxDegreeIndex = Math.Max(maxDegreeIndex, specialGraphDegrees[i].Item1);
                }
                // Initialize the list with default value -1.
                _specialGraphDegrees = Enumerable.Repeat(-1, maxDegreeIndex + 1).ToList();
                // Assign specified degrees to their corresponding levels.
                for (int i = 0; i < specialGraphDegrees.Count; i++)
                {
                    _specialGraphDegrees[specialGraphDegrees[i].Item1] = specialGraphDegrees[i].Item2;
                }
            }
            else
            {
                _specialGraphDegrees = new List<int>();
            }

            // Position the nodes of the inner polygon evenly around a circle.
            UpdateNodesPositions(0, _innerPolygon, 0);
            // Set indices for nodes that have not been updated.
            _pendingNodesStart = 0;
            _pendingNodesEnd = _innerPolygon; // Initially, all inner polygon nodes are "not updated"
        }

        /// <summary>
        /// Copy constructor for TessellationGraph.
        /// Creates a deep copy of the provided TessellationGraph.
        /// </summary>
        /// <param name="other">Another TessellationGraph instance to copy from.</param>
        public TessellationGraph(TessellationGraph other)
        {
            // Copy simple value types.
            _innerPolygon = other._innerPolygon;
            _outerPolygon = other._outerPolygon;
            _defaultGraphDegree = other._defaultGraphDegree;
            _specialGraphDegrees = new List<int>(other._specialGraphDegrees); // Deep copy.
            _radius = other._radius;
            _endLevel = other._endLevel;
            _pendingNodesStart = other._pendingNodesStart;
            _pendingNodesEnd = other._pendingNodesEnd;

            // Deep copy the list of nodes using each Node's copy constructor.
            _nodes = new List<Node>();
            foreach (var node in other._nodes)
            {
                _nodes.Add(new Node(node));
            }

            // Deep copy the polygons list. Each polygon is a list of ints.
            _polygons = new List<List<int>>();
            foreach (var innerList in other._polygons)
            {
                _polygons.Add(new List<int>(innerList));
            }
        }

        // -------------------- Public Methods --------------------

        /// <summary>
        /// Expands the graph by adding new levels (layers).
        /// For each new level, new nodes are created, connected to previous nodes,
        /// and new polygons (faces) are formed.
        /// </summary>
        /// <param name="newLevelCount">How many new levels to add.</param>
        public void UpdateGraph(int newLevelCount)
        {
            // Loop through each new level.
            for (int level = _endLevel; level < _endLevel + newLevelCount; level++)
            {
                // Determine the degree to use for the current level.
                int degreeNow = (level < _specialGraphDegrees.Count && _specialGraphDegrees[level] != -1) ?
                    _specialGraphDegrees[level] : _defaultGraphDegree;

                // Calculate how many new nodes need to be added and update the connection situations.
                int newNodesToBeAdded = CountNewNodesToBeAddedAndUpdateSituations(degreeNow);

                // Create new nodes for the current level.
                for (int i = 0; i < newNodesToBeAdded; i++)
                {
                    _nodes.Add(new Node(level));
                }

                // Perform basic initialization for the new nodes (neighbors list, before/next indices).
                NodesBasicUpdate(newNodesToBeAdded);

                // For each parent node, determine its leftmost and rightmost children.
                SetRightmostLeftmostChildrenInParents(newNodesToBeAdded, degreeNow);

                // For each node in the previous level, assign parents for their new child nodes.
                for (int index = _pendingNodesStart; index < _pendingNodesEnd; index++)
                {
                    SetParentsForNewNodes(index);
                }

                // Update connection metrics (lengths and bridges) for each new child node.
                for (int index = _pendingNodesStart; index < _pendingNodesEnd; index++)
                {
                    SetLenAndBridgeForChildren(index);
                }
                // For new nodes with a special parent marker (-2), form a new polygon.
                for (int pointer = _pendingNodesEnd; pointer < _pendingNodesEnd + newNodesToBeAdded; pointer++)
                {
                    if (_nodes[pointer].Parent == -2)
                        FormPolygon(_nodes[pointer].ParentRight, _nodes[pointer].ParentLeft, pointer);
                }

                // Update the indices for not-yet-updated nodes, these are the nodes in the last layer.
                _pendingNodesStart = _pendingNodesEnd;
                _pendingNodesEnd += newNodesToBeAdded;
                // Update the positions of the newly added nodes.
                UpdateNodesPositions(_pendingNodesStart, _pendingNodesEnd, level);
            }
            // Update the overall level count.
            _endLevel += newLevelCount;
        }

        /// <summary>
        /// Returns a list of polylines representing the tessellated polygons (faces).
        /// </summary>
        /// <returns>List of Polyline objects corresponding to each polygon.</returns>
        public List<Polyline> PolygonsOf()
        {
            var ret = new List<Polyline>();
            // For each polygon, convert node indices to actual points.
            foreach (var polygon in _polygons)
            {
                var points = new List<Point3d>();
                foreach (int index in polygon)
                {
                    points.Add(_nodes[index].Position);
                }
                ret.Add(new Polyline(points));
            }
            return ret;
        }

        /// <summary>
        /// Returns a list of lines (edges) in the graph.
        /// Each line connects two nodes.
        /// </summary>
        /// <returns>List of Line objects representing graph edges.</returns>
        public List<Line> LinesOf()
        {
            var ret = new List<Line>();
            // Loop through each node.
            for (int i = 0; i < _nodes.Count; i++)
            {
                // For each neighbor of the node...
                for (int j = 0; j < _nodes[i].Neighbours.Count; j++)
                {
                    // Only add each edge once.
                    if (_nodes[i].Neighbours[j] > i)
                    {
                        ret.Add(new Line(_nodes[i].Position, _nodes[_nodes[i].Neighbours[j]].Position));
                    }
                }
            }
            return ret;
        }

        /// <summary>
        /// Returns a list of points representing the positions of all nodes in the graph.
        /// </summary>
        /// <returns>List of Point3d objects for each node.</returns>
        public List<Point3d> PointsOf()
        {
            var ret = new List<Point3d>();
            foreach (Node tmp in _nodes)
            {
                ret.Add(tmp.Position);
            }
            return ret;
        }

        /// <summary>
        /// Returns the connectivity information for each polygon.
        /// Each polygon is represented as an array of node indices.
        /// </summary>
        /// <returns>2D array of integers with polygon connectivity.</returns>
        public int[][] PolygonsGraphOf()
        {
            // or: public int[][] GraphOf() => _polygons.Select(innerList => innerList.ToArray()).ToArray();
            int[][] ret = new int[_polygons.Count][];
            for (int i = 0; i < _polygons.Count; i++)
            {
                ret[i] = _polygons[i].ToArray();
            }
            return ret;
        }

        /// <summary>
        /// Returns the connectivity graph for nodes.
        /// Each node's list of neighbor indices is returned.
        /// </summary>
        /// <returns>2D array of integers representing node neighbors.</returns>
        public int[][] GraphOf()
        {
            int[][] ret = new int[_nodes.Count][];
            for (int i = 0; i < _nodes.Count; i++)
            {
                ret[i] = _nodes[i].Neighbours.ToArray();
            }
            return ret;
        }

        /// <summary>
        /// Returns an array representing the layer (level) of each node.
        /// </summary>
        /// <returns>Array of integers where each entry is the node's layer.</returns>
        public int[] LayersOf()
        {
            int[] result = new int[_nodes.Count];
            for (int i = 0; i < _nodes.Count; i++)
            {
                result[i] = _nodes[i].Layer;
            }
            return result;
        }

        // -------------------- Private Helper Methods --------------------

        /// <summary>
        /// Calculates the number of new nodes to be added based on the current level's degree.
        /// Also updates each node's situation regarding its left/right neighbor relationships.
        /// </summary>
        /// <param name="degreeNow">The degree required at the current level.</param>
        /// <returns>The total number of new nodes to add.</returns>
        private int CountNewNodesToBeAddedAndUpdateSituations(int degreeNow)
        {
            int newNodesToBeAdded = 0;
            // Process each node that hasn't been updated yet.
            for (int index = _pendingNodesStart; index < _pendingNodesEnd; index++)
            {
                // Get the previous node index.
                int leftNeighbor = _nodes[index].BeforeIndex;
                // Calculate the current available sides between a node and its left neighbor.
                int sidesWeHave = _nodes[index].LeftLen + _nodes[leftNeighbor].RightLen + _nodes[index].LeftBridge;
                int situation;
                // Determine connection situation based on available sides relative to OuterPolygon.
                if (sidesWeHave + 1 == _outerPolygon)
                {
                    // They need to connect directly.
                    situation = 1;
                }
                else if (sidesWeHave + 2 == _outerPolygon)
                {
                    // They should share a child.
                    situation = 2;
                }
                else
                {
                    // Normal case: no shared child.
                    situation = 3;
                }
                // Update situation values for both the node and its left neighbor.
                _nodes[index].LeftSituation = situation;
                _nodes[leftNeighbor].RightSituation = situation;
                // Increment count of new nodes needed.
                newNodesToBeAdded += degreeNow - _nodes[index].Dimension - (3 - situation); // 3: no overlap, 2: one overlap, 1: two overlaps.
            }
            return newNodesToBeAdded;
        }

        /// <summary>
        /// Performs basic initialization for new nodes added in the update:
        /// sets up their neighbors list and circular 'next' and 'before' indices.
        /// </summary>
        /// <param name="newNodesToBeAdded">Number of new nodes added in the current level.</param>
        private void NodesBasicUpdate(int newNodesToBeAdded)
        {
            // determine the range of new nodes, both inclusive to be more readable
            int firstIndex = _pendingNodesEnd, lastIndex = _pendingNodesEnd + newNodesToBeAdded - 1;
            for (int i = firstIndex; i <= lastIndex; i++)
            {
                // Initialize neighbor list for the new node.
                _nodes[i].Neighbours = new List<int>();
                // Set the 'next' pointer in the circular list.
                _nodes[i].NextIndex = (i != lastIndex) ? i + 1 : firstIndex;
                // Set the 'before' pointer.
                _nodes[i].BeforeIndex = (i != firstIndex)? i - 1: lastIndex;
            }
        }


        /// <summary>
        /// Determines and assigns the leftmost and rightmost children for each node (parent)
        /// based on how many new child nodes they should have.
        /// </summary>
        /// <param name="NewNodesToBeAdded">Number of new nodes in the current level.</param>
        /// <param name="DegreeNow">The required degree for the current level.</param>
        private void SetRightmostLeftmostChildrenInParents(int NewNodesToBeAdded, int DegreeNow)
        {
            int pointer = _pendingNodesEnd; // Start pointer at the beginning of new nodes.
            // Iterate over each parent node in the previous level.
            for (int index = _pendingNodesStart; index < _pendingNodesEnd; index++)
            {
                // Calculate how many children this node should have.
                int childrenThisNodeShouldHave = DegreeNow - _nodes[index].Dimension;
                // If we have a connection to the left or right node, no child corresponding to them should be added
                if (_nodes[index].LeftSituation == 1) childrenThisNodeShouldHave--;
                if (_nodes[index].RightSituation == 1) childrenThisNodeShouldHave--;
                // set the LeftMostChild
                if (index == _pendingNodesStart && _nodes[index].LeftSituation == 2)
                {
                    // Special case for the first node if has a shared child with last node (circling over the array)
                    _nodes[index].LeftmostChild = _pendingNodesEnd + NewNodesToBeAdded - 1;
                    childrenThisNodeShouldHave--; // Adjust due to special case.
                }
                else
                {
                    // Normal situation
                    _nodes[index].LeftmostChild = pointer;
                }
                // Advance the pointer to mark the rightmost child.
                for (int i = 0; i < childrenThisNodeShouldHave - 1; i++) pointer = _nodes[pointer].NextIndex;
                _nodes[index].RightmostChild = pointer;

                // Adjust pointer to the next child if the right side is not a shared child.
                if (_nodes[index].RightSituation != 2)
                {
                    pointer = _nodes[pointer].NextIndex;
                }
            }
        }

        /// <summary>
        /// Assigns parent relationships for new nodes based on their position between leftmost and rightmost children.
        /// Also handles special cases where nodes share a parent.
        /// </summary>
        /// <param name="index">The index of the parent node for which to assign children.</param>
        private void SetParentsForNewNodes(int index)
        {
            int pointer = _nodes[index].LeftmostChild;
            // Loop through all new children of the current parent node.
            while (true)
            {
                // By default, assign the parent.
                _nodes[pointer].Parent = index;
                // Handle special situation for left shared child.
                if (pointer == _nodes[index].LeftmostChild && _nodes[index].LeftSituation == 2)
                {
                    _nodes[pointer].Parent = -2;  // Marker indicating shared parent.
                    _nodes[pointer].ParentRight = index;
                }
                // Handle special situation for right shared child.
                if (pointer == _nodes[index].RightmostChild && _nodes[index].RightSituation == 2)
                {
                    _nodes[pointer].Parent = -2;  // Marker indicating shared parent.
                    _nodes[pointer].ParentLeft = index;
                }

                // Break once the last child is reached.
                if (pointer == _nodes[index].RightmostChild)
                {
                    break;
                }
                else
                {
                    pointer = _nodes[pointer].NextIndex;
                }
            }
        }

        /// <summary>
        /// Updates the connection metrics (length and bridge) for each child node,
        /// and connects the child nodes back to their parent.
        /// </summary>
        /// <param name="index">The index of the parent node to process.</param>
        private void SetLenAndBridgeForChildren(int index)
        {
            int pointer = _nodes[index].LeftmostChild;
            // Special case: if left situation requires direct connection.
            if (_nodes[index].LeftSituation == 1)
            {
                Connect(index, _nodes[index].BeforeIndex);
                // Form a polygon between the current node and its neighbor.
                FormPolygon(index, _nodes[index].BeforeIndex, -1);
            }
            // Loop through all child nodes.
            while (true)
            {
                // Connect the child node with its parent.
                Connect(pointer, index);
                // Set default metrics for the child node.
                _nodes[pointer].LeftLen = 1;
                _nodes[pointer].RightLen = 1;
                _nodes[pointer].LeftBridge = 0;
                _nodes[pointer].RightBridge = 1;
                // Adjust metrics based on the left-most child situation.
                if (pointer == _nodes[index].LeftmostChild)
                {
                    if (_nodes[index].LeftSituation == 1)
                    {
                        _nodes[pointer].LeftLen = 1;
                        _nodes[pointer].LeftBridge = 1;
                    }
                    else if (_nodes[index].LeftSituation == 2)
                    {
                        // Note: Degree should be greater than 3 to avoid bugs.
                        _nodes[pointer].LeftLen = 1;
                        _nodes[pointer].LeftBridge = 0;
                    }
                    else if (_nodes[index].LeftSituation == 3)
                    {
                        _nodes[pointer].LeftLen = _nodes[index].LeftLen + 1;
                        _nodes[pointer].LeftBridge = _nodes[index].LeftBridge;
                    }
                }
                // Adjust metrics for the right-most child situation.
                if (pointer == _nodes[index].RightmostChild)
                {
                    if (_nodes[index].RightSituation == 1)
                    {
                        _nodes[pointer].RightLen = 1;
                        _nodes[pointer].RightBridge = 1;
                    }
                    else if (_nodes[index].RightSituation == 2)
                    {
                        // Note: Degree should be greater than 3 to avoid bugs.
                        _nodes[pointer].RightLen = 1;
                        _nodes[pointer].RightBridge = 0;
                    }
                    else if (_nodes[index].RightSituation == 3)
                    {
                        _nodes[pointer].RightLen = _nodes[index].RightLen + 1;
                        _nodes[pointer].RightBridge = _nodes[index].RightBridge;
                    }
                }
                // Break if the last child has been processed. And if not, go to the next node
                if (pointer == _nodes[index].RightmostChild)
                {
                    break;
                }
                else
                {
                    pointer = _nodes[pointer].NextIndex;
                }
            }
        }

        /// <summary>
        /// Constructs a new polygon (face) by collecting node indices from two branches.
        /// It follows the left and right lengths from the given nodes and combines them.
        /// </summary>
        /// <param name="right">Index of the right reference node.</param>
        /// <param name="left">Index of the left reference node.</param>
        /// <param name="root">
        /// If not -1, the root node index is appended to the polygon.
        /// </param>
        private void FormPolygon(int right, int left, int root)
        {
            var rightBranch = new List<int>(); // List for the right branch.
            var leftBranch = new List<int>();  // List for the left branch.
            int rightPointer = right, leftPointer = left;
            // Traverse the right branch for the number of steps given by LeftLen.
            for (int i = 0; i < _nodes[right].LeftLen; i++)
            {
                rightBranch.Add(rightPointer);
                // If the node is marked as shared, follow the special pointer.
                if (_nodes[rightPointer].Parent == -2)
                    rightPointer = _nodes[rightPointer].ParentLeft;
                else
                    rightPointer = _nodes[rightPointer].Parent;
            }
            // Traverse the left branch for the number of steps given by RightLen.
            for (int i = 0; i < _nodes[left].RightLen; i++)
            {
                leftBranch.Add(leftPointer);
                if (_nodes[leftPointer].Parent == -2)
                    leftPointer = _nodes[leftPointer].ParentRight;
                else
                    leftPointer = _nodes[leftPointer].Parent;
            }
            // Reverse the left branch so that the polygon points follow a continuous order.
            leftBranch.Reverse();
            // Combine the branches.
            rightBranch.Add(rightPointer);
            if (rightPointer != leftPointer)
                rightBranch.Add(leftPointer);
            rightBranch.AddRange(leftBranch);
            // Append the root node if provided.
            if (root != -1)
                rightBranch.Add(root);
            // Add the constructed polygon to the list.
            _polygons.Add(rightBranch);
        }

        /// <summary>
        /// Updates the positions of nodes in a given range.
        /// Nodes are placed evenly along a circle based on the current level.
        /// </summary>
        /// <param name="start">Starting index of nodes to update.</param>
        /// <param name="end">Ending index (exclusive) of nodes to update.</param>
        /// <param name="level">The current level (layer) for positioning.</param>
        private void UpdateNodesPositions(int start, int end, int level)
        {
            for (int i = start; i < end; i++)
            {
                // Calculate angle (theta) for even spacing.
                double teta = (Math.PI * 2 / (end - start)) * (i - start);
                // Compute the position on a circle with a radius that scales with the level.
                _nodes[i].Position = new Point3d(Math.Cos(teta) * (level + 1) * _radius,
                                            Math.Sin(teta) * (level + 1) * _radius,
                                            0);
            }
        }

        /// <summary>
        /// Connects two nodes by adding each node to the other's neighbors list,
        /// and increments their connection (dimension) counters.
        /// </summary>
        /// <param name="node1">Index of the first node.</param>
        /// <param name="node2">Index of the second node.</param>
        private void Connect(int node1, int node2)
        {
            _nodes[node1].Neighbours.Add(node2);
            _nodes[node2].Neighbours.Add(node1);
            _nodes[node1].Dimension++;
            _nodes[node2].Dimension++;
        }

        // -------------------- Nested Node Class --------------------

        /// <summary>
        /// The Node class represents a single vertex in the tessellation graph.
        /// It stores its position, connectivity data, and other metadata used
        /// during the tessellation process.
        /// </summary>
        class Node
        {
            // -------------------- Public Fields --------------------

            /// <summary>
            /// The position of the node.
            /// </summary>
            public Point3d Position;

            /// <summary>
            /// Index of the next node in the circular order.
            /// </summary>
            public int NextIndex;

            /// <summary>
            /// Index of the previous node in the circular order.
            /// </summary>
            public int BeforeIndex;

            /// <summary>
            /// The current number of connections (neighbors) the node has.
            /// </summary>
            public int Dimension;

            /// <summary>
            /// Length (or count) of connected edges on the left side.
            /// </summary>
            public int LeftLen;

            /// <summary>
            /// Length (or count) of connected edges on the right side.
            /// </summary>
            public int RightLen;

            /// <summary>
            /// A 0 or 1 metric for the left bridge connection.
            /// if it is one, means that the top node on the branch has a 
            /// connection to a node in the same level (a bridge)
            /// </summary>
            public int LeftBridge;

            /// <summary>
            /// A 0 or 1 metric for the right bridge connection.
            /// if it is one, means that the top node on the branch has a 
            /// connection to a node in the same level (a bridge)
            /// </summary>
            public int RightBridge;

            /// <summary>
            /// The layer (level) of this node in the graph.
            /// Not useful in shaping the graph, but a useful output for user.
            /// </summary>
            public int Layer;

            /// <summary>
            /// The parent node index for hierarchical relationships.
            /// If it is -2, means that the node has two parents, right and left
            /// </summary>
            public int Parent;

            /// <summary>
            /// For shared nodes, the right parent's index.
            /// </summary>
            public int ParentLeft;

            /// <summary>
            /// For shared nodes, the left parent's index.
            /// </summary>
            public int ParentRight;

            /// <summary>
            /// The situation for the left neighbor connection.
            /// 1: direct connect (bridge), 2: share a child, 3: normal.
            /// </summary>
            public int LeftSituation;

            /// <summary>
            /// The situation for the right neighbor connection.
            /// 1: direct connect (bridge), 2: share a child, 3: normal.
            /// </summary>
            public int RightSituation;

            /// <summary>
            /// Index of the leftmost child in the hierarchical connection.
            /// Also counting shared children
            /// </summary>
            public int LeftmostChild;

            /// <summary>
            /// Index of the rightmost child in the hierarchical connection.
            /// Also counting shared children
            /// </summary>
            public int RightmostChild;

            /// <summary>
            /// List of neighbor node indices.
            /// </summary>
            public List<int> Neighbours;

            // -------------------- Constructors --------------------

            /// <summary>
            /// Constructor for nodes in the first (inner) layer.
            /// Sets initial values and neighbor indices based on the inner polygon.
            /// </summary>
            /// <param name="before">Index of the previous node in the inner polygon.</param>
            /// <param name="next">Index of the next node in the inner polygon.</param>
            public Node(int before, int next)
            {
                Position = new Point3d(0, 0, 0);   // Initialize position.
                NextIndex = next;
                BeforeIndex = before;
                Neighbours = new List<int>(); // Initialize empty neighbors list.
                Dimension = 0;                // No connections yet.
                Layer = 0;                    // First layer.
                LeftLen = RightLen = 0;
                LeftBridge = RightBridge = 1; // Default bridge values.
                LeftSituation = RightSituation = -1;      // Uninitialized situation.
                LeftmostChild = RightmostChild = -1;
                Parent = ParentLeft = ParentRight = -1;
            }

            /// <summary>
            /// Constructor for nodes in subsequent layers.
            /// </summary>
            /// <param name="layer">The layer (level) this node belongs to.</param>
            public Node(int layer)
            {
                Position = new Point3d(0, 0, 0);   // Initialize position.
                NextIndex = 0;
                BeforeIndex = 0;
                Neighbours = new List<int>(); // Initialize neighbors list.
                Layer = layer;
                LeftLen = RightLen = 0;
                LeftBridge = RightBridge = 0;
                Dimension = 0;
                LeftSituation = RightSituation = 0;
                LeftmostChild = RightmostChild = -1;
                Parent = ParentLeft = ParentRight = -1;
            }

            /// <summary>
            /// Copy constructor for Node.
            /// Creates a deep copy of the given node.
            /// </summary>
            /// <param name="other">Another Node instance to copy from.</param>
            public Node(Node other)
            {
                // Copy all value-type fields.
                Position = other.Position; // Point3d is a struct.
                NextIndex = other.NextIndex;
                BeforeIndex = other.BeforeIndex;
                Dimension = other.Dimension;
                LeftLen = other.LeftLen;
                RightLen = other.RightLen;
                LeftBridge = other.LeftBridge;
                RightBridge = other.RightBridge;
                Layer = other.Layer;
                Parent = other.Parent;
                ParentLeft = other.ParentLeft;
                ParentRight = other.ParentRight;
                LeftSituation = other.LeftSituation;
                RightSituation = other.RightSituation;
                LeftmostChild = other.LeftmostChild;
                RightmostChild = other.RightmostChild;
                // Deep copy the list of neighbors.
                Neighbours = new List<int>(other.Neighbours);
            }
        } // End of Node class.
    } // End of TessellationGraph class.
}