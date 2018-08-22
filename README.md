# graph-signal-processing 
GRAPH signal processing (GSP) is an emerging field that studies signals supported on irregular domains. It extends traditional signal processing techniques to more intricate data structures, finding applications in sensor networks, image processing, and clustering, to name a few. This project applies graph signal processing to the semi-supervised field of handwritten numeral recognition to achieve handwritten numeral classification.  
notions：  
1.graph: Graphs are generic data representation forms which are useful for describing the geometric structures of data domains in numerous applications, including social, energy, transportation, sensor and neuronal networks. (Non-rule discrete domains)//
2.A graph signal: A graph signal is a real valued function defined on each node of the graph//
3.sampling theorem of graph signal//
What? ——Sampling theory of graph signals similarly deals with the problem of recovering a signal from its samples on a subset of nodes of the graph.//
Why? ——The sampling theory for graph signals aims to extend the traditional Nyquist-Shannon sampling theory by allowing us to identify the class of graph signals that can be reconstructed from their values on a subset of vertices.//
Apply://
1.In this experiment, we used our proposed active semisupervised learning algorithm to perform a classification task on the USPS handwritten digits dataset2. This dataset consists of 1100 16*16 pixel images for each of the digits 0
to 9. We used 100 randomly selected samples for each digit class to create one instance of our dataset. Thus each instance consists of 1000 feature vectors of dimension 256. The graph is constructed using Gaussian kernel weights, where xi is the 256-dimensional
feature vector composed of pixel intensity values for each image. The parameter  is chosen to be 1=3-rd of the average distance to the K-th nearest neighbor for all datapoints. This heuristic has been suggested in. Additionally, the graph is sparsied approximately by restricting
the connectivity of each datapoint to its K nearest neighbors, i.e., an edge between nodes i and j is removed unless node i is among the K-nearest neighbors of node j or vice-versa. This results in a symmetric adjacency matrix for the graph. Using the graph constructed, we select
the points to label and report prediction error after reconstruction using our semi-supervised learning algorithm. We repeat the classication over 10 such instances of the dataset
and report the average classification error. The results are illustrated in Figure (5a). We observe that our proposed method outperforms the others. A notable feature of our method is that we show very good classification results even for very few labeled samples. This is due to our inherent
criterion for active learning that tries to select those points that maximize the recoverable dimensions of the underlying data manifold.
