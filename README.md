# RF-clustering

### Step 1: Create synthetic label
..* Addcl1: ‘synthetic-labeled’ data are added by randomly sampling from the product of empirical marginal distributions of the variables.
The Addcl1 RF dissimilarity weighs the contribution of each variable on the dissimilarity according to how dependent it is on the other variables.

..* Addcl2: ‘synthetic-labeled’ data are added by randomly sampling from the hyper rectangle that contains observed data. (Suitable when MDS plot leads to two distinct point clouds corresponding to the values of the binary data.)


### Step 2: Perform Random Forest Predictor
..* Number of forests: 100
..* Number of trees: 4000
..* The idea is to use the similarity matrix generated from a RF predictor that distinguishes observed from “synthetic” data.
The RF dissimilarity easily deals with a large number of variables due to its intrinsic variable selection.
..* Here we use Addcl1 proximity to calculate by 𝑑𝑖𝑠𝑠𝑖𝑚𝑖𝑙𝑎𝑟𝑖𝑡𝑦=√(1−𝑝𝑟𝑜𝑥𝑖𝑚𝑖𝑡𝑦) (Proximity measures among the input based on the frequency that pairs of data points are in the same terminal nodes)


### Step 3: Approximating the RF dissimilarity
..* When dealing with quantitative variables, one can sometimes find a Euclidean distance-based approximation of the Addcl1 RF dissimilarity if each variable is equally important for distinguishing observed from synthetic observations.
..* RF dissimilarity depends only on variable ranks since the underlying tree node splitting criterion (Gini index) considers only variable ranks
..* Using the resulting variables in a Euclidean distance


### Classical MDS for Addcl1 dissimilarity
..* The RF dissimilarity can be used as input of MDS, which yields a set of points in an Euclidean space such that the Euclidean distances between these points are approximately equal to the dissimilarities.

..* Multidimensional scaling (MDS) algorithms start with a matrix of item-item distances and then assign coordinates for each item in a low-dimensional space to represent the distances graphically. Unlike other ordination methods, MDS makes few assumptions about the nature of the data. For example, principal components analysis assumes linear relationships and reciprocal averaging assumes modal relationships. MDS makes neither of these assumptions, so is well suited for a wide variety of data.

..* Usually, we use classical MDS for the Addcl1 dissimilarity to avoid the detection of spurious patterns. There is empirical evidence that the Addcl1 RF dissimilarity can be superior to standard distance distance measures in several applications.


### Next steps
1. Reduce Spearman rank correlation threshold in Step 1
2. Examine imaging features of clusters
..* Imaging features used in analysis
..* Imaging features not used in analysis
3. Compare clusters based on their
..* Clinical characteristics
..* RNA sequencing data (when available)
