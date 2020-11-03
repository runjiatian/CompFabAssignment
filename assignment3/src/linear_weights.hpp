#include <Eigen/Core>
#include <limits>

// TODO: HW3
// Assignment 3, Part 2.1. 
// Fill the values inside linear weight matrix W of shape (n, m)
// V is a vertices matrix in shape (n, 3), each row is the position of a vertex in mesh
// C is a control points matrix in shape (m, 3), each row is the position of a control point
void ComputeLinearSkinningWeights(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXd& C, 
	Eigen::MatrixXd& W) 
{
	W.resize(V.rows(), C.rows());

	/* Implement your code here. */
    // Loop through the vertices matrix
    for (int i = 0; i < V.rows(); ++i) {
        double sum = 0;
        // Loop through the handles matrix
        for (int j = 0; j < C.rows(); ++j) {
            double dist = (V.row(i) - C.row(j)).norm();
            W(i, j) = dist;
            sum += dist;
        }

        // Normalize with sum
        for (int j = 0; j < C.rows(); ++j) {
            W(i, j) = 1/W(i, j);
        }
    }
}