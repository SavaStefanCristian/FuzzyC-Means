# Fuzzy C-Means Clustering in C++

## Overview
This repository contains an implementation of the **Fuzzy C-Means (FCM) clustering algorithm** in C++. The program clusters randomly generated points in an n-dimensional Euclidean space based on their similarity.

## Features
- Implements **Fuzzy C-Means (FCM) algorithm** for clustering.
- Uses **Euclidean distance** for computing cluster memberships.
- Supports **custom fuzziness parameters** to control clustering behavior.
- Dynamically generates **random data points** for testing.
- Outputs **centroids and membership matrices** to a file.

## Algorithm Complexity
The upper asymptotic complexity of this FCM algorithm is approximately:

```
O(K * N * D * C)
```
Where:
- **K** = Number of iterations until convergence
- **N** = Number of data points
- **D** = Number of dimensions
- **C** = Number of clusters


#### Example Input
```
Enter the number of dimensions (n): 2
Enter the number of randomly generated points: 10
Enter the number of clusters: 3
Enter the fuzziness parameter: 2.0
```
