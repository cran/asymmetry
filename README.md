## Asymmetry

This R package provides a comprehensive suite of tools for the analysis and visualization of asymmetric data using Multidimensional Scaling (MDS) and related decomposition techniques. Asymmetry is a common challenge in social sciences, marketing, and sensory research, where the relationship from object $i$ to $j$ does not necessarily equal the relationship from $j$ to $i$.
The package implements several specialized models to capture these directional differences:

ASYMSCAL: A weighted MDS model that accounts for asymmetry through differential weighting of dimensions, estimated via the SMACOF (Scaling by Majorizing a Complicated Function) framework.

MDS with Unique Dimensions: An extension that identifies specific dimensions or "uniqueness" values unique to individual rows or columns within a distance matrix.

Slide-vector Model: A method that models asymmetry by adding a directional "slide" vector to a symmetric configuration, representing a shared bias or trend in the data.

Skew-symmetric Decomposition: Tools to decompose a square matrix into its symmetric and skew-symmetric components ($S$ and $A$), including singular value decomposition (SVD) of the skew-symmetric part for low-dimensional representation.
Beyond model fitting, the package includes robust diagnostic and visualization functions. It offers heatmap visualizations for dominance relations, as well as specialized plot methods for configurations, slide-vectors, and skew-symmetric planes.
