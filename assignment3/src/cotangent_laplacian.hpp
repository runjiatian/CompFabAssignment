#pragma once

#include <Eigen/Sparse>
#include <igl/edge_lengths.h>
#include <igl/face_areas.h>
#include <igl/dihedral_angles.h>
#include <igl/volume.h>
#include "dihedral_sine.hpp"

template <typename T>
using Vector3 = Eigen::Matrix<T, 3, 1>;

// TODO: HW3
// Assignment 3, Part 3.2.
/* Implement your code here. */
// Implement the function to compute the cotangent laplacian matrix L 
// V is the vertex matrix of shape (n, 3), each row is the position of a vertex in the mesh
// F is the element index matrix of shape (m, 4), each row is the vertex indices of a tetrahedron
// L is the output cotangent laplacian matrix of shape (n, n), and it's a sparse matrix.
// Hints:
	// 1. For each tetrahedron, loop over each of its edge,
	//    consider which part of the L matrix this edge in this tetrahedron contributes to
	// 2. compute the cos and sin of the dihedral angle by the law of diehedral angles http://mathworld.wolfram.com/Tetrahedron.html
	//	  specifically, compute the sin and cos of dihedral angles from the edge lengths, face areas and tet volume
	// 3. build the triplets <row, col, value> in IJV
void cotangent_laplacian(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F, 
	Eigen::SparseMatrix<double>& L) 
{
	L.resize(V.rows(), V.rows());

	std::vector<Eigen::Triplet<double> > IJV;
	IJV.clear();

	/* Implement your code here. */
	// Initialize weight to be zero
    Eigen::MatrixXd Length_o;
    Eigen::MatrixXd Area;
    Eigen::MatrixXd Volume;

//    //[3 0],[3 1],[3 2],[1 2],[2 0],[0 1]
//    igl::edge_lengths(V, F, Length_o);
//
//    //by 4 list of face areas corresponding to faces opposite vertices
//    //0,1,2,3
//    igl::face_areas(V, F, Area);
//    igl::volume(V, F, Volume);

//    Eigen::MatrixXd Length = Length_o;
//
//    Length.col(0) = Length_o.col(3);
//    Length.col(1) = Length_o.col(4);
//    Length.col(2) = Length_o.col(5);
//    Length.col(3) = Length_o.col(0);
//    Length.col(4) = Length_o.col(1);
//    Length.col(5) = Length_o.col(2);
//
//    //https://github.com/libigl/libigl/blob/master/include/igl/face_areas.h
//    Eigen::Matrix<double, Eigen::Dynamic, 6> Thr;

    //Thr computes sine value for edge 12, 20,01,30,31,32
    //dihedral_sine(Volume, Area, Length, Thr);


    //Initialize to Zero
	// Loop through vertex
	// Loop through edges
	    //For source and destination, store the value of cotangent
        //For destination and destination, store the value of cotangent
        //Increment Lii for i != j
        //
    for (int r = 0; r < F.rows(); ++r) {
        int v1 = F(r, 0);
        int v2 = F(r, 1);
        int v3 = F(r, 2);
        int v4 = F(r, 3);

        Eigen::MatrixXd w(4,4);

        Vector3<double> V1 = V.row(v1);
        Vector3<double> V2 = V.row(v2);
        Vector3<double> V3 = V.row(v3);
        Vector3<double> V4 = V.row(v4);

        double l12 = (V1 - V2).norm();
        double l13 = (V1 - V3).norm();
        double l14 = (V1 - V4).norm();
        double l23 = (V3 - V2).norm();
        double l24 = (V4 - V2).norm();
        double l34 = (V4 - V3).norm();

        // Computer Innter Pointing Triangle Normal
        Vector3<double> n134 = ((V1 - V3).cross(V4 - V3))/(((V1 - V3).cross(V4 - V3)).norm());
        Vector3<double> n124 = ((V2 - V1).cross(V4 - V1))/(((V2 - V1).cross(V4 - V1)).norm());
        Vector3<double> n234 = ((V4 - V3).cross(V2 - V3))/(((V4 - V3).cross(V2 - V3)).norm());
        Vector3<double> n123 = ((V3 - V1).cross(V2 - V1))/(((V3 - V1).cross(V2 - V1)).norm());

        // Compute Angle
        // Edge 12
        w(0,1) = (-n134.dot(n234))/sqrt(1 - n134.dot(n234)*n134.dot(n234));
        // Edge 13
        w(0,2) = (-n124.dot(n234))/sqrt(1 - n124.dot(n234)*n124.dot(n234));
        // Edge 14
        w(0,3) = (-n123.dot(n234))/sqrt(1 - n123.dot(n234)*n123.dot(n234));
        //Edge 23
        w(1,2)  = (-n134.dot(n124))/sqrt(1 - n134.dot(n124)*n134.dot(n124));
        // Edge 24
        w(1,3)  = (-n123.dot(n134))/sqrt(1 - n123.dot(n134)*n123.dot(n134));
        // Edge 34
        w(2,3) = (-n123.dot(n124))/sqrt(1 - n123.dot(n124)*n123.dot(n124));

        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j){
                int vi = F(r, i);
                int vj = F(r, j);
                int a = std::min(i, j);
                int b = std::max(i, j);
                IJV.push_back(Eigen::Triplet<double>(vi, vj, w(a,b)));
                IJV.push_back(Eigen::Triplet<double>(vj, vi, w(a,b)));
                IJV.push_back(Eigen::Triplet<double>(vi, vi, -w(a,b)));
                IJV.push_back(Eigen::Triplet<double>(vj, vj, -w(a,b)));
            }
    }

	// Set From Triplets Sums all Triplets with the same indices
	L.setFromTriplets(IJV.begin(), IJV.end());
}