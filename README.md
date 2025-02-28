# Hyperbolic Tessellations
### Using Physics-Based Relaxation in Form-Finding of Non-Euclidean Surfaces

**Personal Project: February 2025**

## Abstract

In this project, I build on my GPU-accelerated truss simulation framework to explore physics-based form-finding of non-Euclidean geometries. By iteratively constructing a planar tessellation in which each vertex connects to a fixed number of neighbors and each polygon maintains a fixed number of edges, I form a system of mass-spring-dampers where diagonals and perimeter beams preserve shape integrity.

Through low damping and dynamic relaxation, these tessellations transition from a flat layout to stable hyperbolic forms, with each polygon remaining nearly planar and edges retaining equal lengths. Leveraging GPU parallelization significantly boosts computation speed, making large-scale form-finding feasible. Additionally, by constraining the outer nodes along a chosen boundary curve, this method can be applied to tent-like architectural designs and other complex surfaces.

## Project Description

### TessellationGraph: Making Initial Planar Graph

Before introducing physics-based relaxation, a planar tessellated graph is generated where each node (vertex) has a fixed number of neighbors (its “degree”) and every face (polygon) has a fixed number of edges. The process begins with an inner polygon—a simple loop of nodes arranged on a circle—that forms the base polygon by connecting adjacent nodes.

Subsequent rings of nodes are then iteratively added around this initial layer. At each new level, the algorithm calculates the exact number of new nodes required between existing ones so that both new and old nodes achieve the correct degree. New nodes are positioned on an expanding circle and connected to the previous layer using a hierarchical parent–child system. Special rules handle edge cases such as shared children or connections between nodes on the same level. In doing so, the algorithm also tracks polygon formation by traversing connected nodes, ensuring that each face maintains the required edge count.

The algorithm supports outer polygons with a fixed number of edges (with the inner polygon possibly differing) and allows for variable node degrees at different levels to introduce variety and complexity into the tessellation.

### Relaxation Phase: From Graph to Form

Once the tessellation graph is generated, the system is reconstructed as a mass-spring-damper model using a custom Truss library. In this phase, nodes move with low damping, allowing the structure to settle into a stable configuration. Key steps include:

- **Diagonals:** Added within each polygon face to help preserve its shape during relaxation.
- **Perimeter Filtering:** Non-polygon nodes and edges are removed, ensuring that only structurally relevant components remain.
- **Helping Beams:** Inserted between neighboring nodes to prevent the structure from folding or collapsing.
- **Diminishing Factors:** Applied on outer layers to subtly reduce their size and stiffness, maintaining aesthetic balance and preventing the outer geometry from overpowering the inner structure.

Together, these elements guide the flat tessellation into a smooth, hyperbolic form that retains nearly equal edge lengths and planar faces.

### Relaxations: Examples

The final relaxed forms demonstrate the power of this approach to achieve complex geometries that would be difficult to obtain through traditional methods. While individual polygons may bend slightly due to compromises in the relaxation settings, overall deformations are minimal. In some cases, polygons from different branches may overlap; these issues can be addressed during subsequent design stages.

This technique shows potential for a wide range of applications—from sculptural elements like lamps to large-scale architectural installations such as expo pavilions.

## Running the Project

All necessary DLLs are provided in the repository. You can run the included examples directly in Rhino/Grasshopper:

- **Rhino/Grasshopper Integration:**  
  Load the provided DLLs into your Rhino/Grasshopper environment.
- **Examples:**  
  Open the example files to see the tessellation graph generation and the subsequent relaxation process in action.


Feel free to explore the code, run the examples, and experiment with the parameters. Contributions and suggestions are welcome!
