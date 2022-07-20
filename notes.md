# Neural networks to learn Boolean (like) networks from multiomics data

## Strategy 1:
Single model per node, no signal integration
Can learn models with good accuracy, independent of in-degree

## Strategy 2:
Integrated model for network - one model for entire network
Aware that GNN will give best results but want to try manipulating features to reuse above model
Input - ragged tensor, only input features from PKN
Target - single node
Multiply input features by the adjacency matrix?

## Strategy 3:
Use tensorflow.gnn package to learn node/edge features and built-in simulation
