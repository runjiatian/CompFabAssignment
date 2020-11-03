#pragma once

#include "read_stl.hpp"
#include "BasicGeometry.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <unordered_set>
#include <fstream>
#include <ctime>
#include <map>
#include <set>
#include "hexahedral_mesh.hpp"

namespace mesh {

	template<typename T>
	class Voxelizer {
	public:
		Voxelizer(const std::string& stl_file_name, const T dx)
			: _dx(dx) {
			// Randomness.
			srand(static_cast<unsigned>(time(0)));
			// Load triangles from the stl file.
			std::vector<Vector3<T>> normals;
			if (!ReadSTL(stl_file_name, _triangles, normals)) {
				std::cout << "ERROR: cannot read " << stl_file_name << std::endl;
				return;
			}
			// Compute the bounding box of _triangle and save the results into _pmin.
			_pmin = _triangles[0][0];
			Vector3<T> pmax = _triangles[0][0];
			for (const auto& triangle : _triangles)
				for (const auto& v : triangle) {
					_pmin = _pmin.cwiseMin(v);
					pmax = pmax.cwiseMax(v);
				}
			for (int i = 0; i < 3; ++i) {
				_pmin[i] -= _dx;
				pmax[i] += _dx;
			}
			// Compute the number of voxels along each direction.
			for (int i = 0; i < 3; ++i)
				_nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;
			// Initialize the voxel array.
			_voxels = std::vector<std::vector<std::vector<bool>>>(_nvoxel.x(),
				std::vector<std::vector<bool>>(_nvoxel.y(),
					std::vector<bool>(_nvoxel.z(), false)));
		}

		const Vector3<T> pmin() const { return _pmin; }
		const T dx() const { return _dx; }
		const Vector3<T> pmax() const { return _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx; }
		const Vector3<int> voxel_num() const { return _nvoxel; }

		// TODO: HW2
		// part 1.2.
		// Fill the _voxels array with the correct flag.
		void BasicVoxelization() {
			const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2]; // number of voxels in each dimension
			for (int i = 0; i < nx; ++i)
				for (int j = 0; j < ny; ++j)
					for (int k = 0; k < nz; ++k)
						_voxels[i][j][k] = false;

            std::vector<T> _hits;
            std::vector<geometry::Triangle<T>> AllTriangles;
            for(int m = 0; m < _triangles.size(); ++m)
            {
                Vector3<T> v0 = _triangles[m][0];
                Vector3<T> v1 = _triangles[m][1];
                Vector3<T> v2 = _triangles[m][2];

                // Step 2: shoot rays.
                geometry::Triangle<T> CurTri(v0, v1, v2);
                AllTriangles.push_back(CurTri);
            }

            // Loop through grirds on XY plane, shoot along z axis
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                {
                    // Step 1: Clear _hits.
                    _hits.clear();

                    // Origin Point of ray
                    T basex = _pmin[0] + (i + 0.5) * _dx;
                    T basey = _pmin[1] + (j + 0.5) * _dx;
                    Vector3<T> base = Vector3<T>(basex, basey, 0);

                    // one way to deal with corner case, everytime have a small pertubation
                    // and counts whether you hit multiple times on a location
                    // e.g. for some hits, more than 1/2 of hitting, then marks the number
                    Vector3<T> dir = Vector3<T>(0,0,1);

                    for(int m = 0; m < AllTriangles.size(); ++m)
                    {
                        // Step 2: shoot rays.
                        geometry::Triangle<T> Curtriangle = AllTriangles[m];
                        Vector3<T> v0 = Curtriangle.vertices(0);
                        Vector3<T> v1 = Curtriangle.vertices(1);
                        Vector3<T> v2 = Curtriangle.vertices(2);

                        Vector3<T> pmin = v2.cwiseMin(v0.cwiseMin(v1));
                        Vector3<T> pmax = v2.cwiseMax(v0.cwiseMax(v1));

                        if((pmin[0] > basex && pmin[1] > basey) || (pmax[0] < basex && pmax[1] < basey)) continue;
                        else{
                            const T t = Curtriangle.IntersectRay(base, dir);
                            if (t == -1) continue;
                            _hits.push_back(base[2] + t * dir[2]);
                        }
                    }

                    if(_hits.empty()) continue;

//                    std::sort(_hits.begin(), _hits.end(),
//                              [](const Vector3<T>A, const Vector3<T>B) ->bool {
//                        if(A[2] - B[2] > 1e-7) return true;
//                        return false;
//                    });

                    std::sort(_hits.begin(),_hits.end(), [](const T A, const T B) -> bool
                    {if (B - A > 1e-7) return true; return false;});

                    // Step 3: fill the _voxels array

                    /*Simpler way: The in-mesh point and out-mesh point, are actually two pairs
                     * Every time increment by 2, get durrent and next index, process enter and exit point,
                     * and loop through only these voxels*/

                    for (int idxIn = 0; idxIn < _hits.size(); idxIn += 2)
                    {
                        int idxOut = idxIn + 1;

                        T zmin = _hits[idxIn];
                        T zmax = _hits[idxOut];

                        int nzmin = static_cast<int> ((zmin - _pmin[2] + 0.5*_dx)/_dx);
                        int nzmax = static_cast<int> ((zmax - _pmin[2]- 0.5*_dx)/_dx);

                        for(int n = nzmin; n < nzmax; ++n)
                        {
                            _voxels[i][j][n] = true;
                        }
                    }
                }

		}

		// TODO: HW2
		// part 2.1.
		// Fill the _voxels array with the correct flag.
		void AdvancedVoxelization() {
			const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2]; // number of voxels in each dimension
            // Step 1: Clear _voxels
            std::vector<std::vector<std::vector<T>>> _hits =
                    std::vector<std::vector<std::vector<T>>>(nx,
                            std::vector<std::vector<T>>(ny,
                                         std::vector<T>(1, -10000000)));

			for (int i = 0; i < nx; ++i)
				for (int j = 0; j < ny; ++j)
					for (int k = 0; k < nz; ++k)
						_voxels[i][j][k] = false;

            for(int m = 0; m < _triangles.size(); ++m)
            {
                Vector3<T> v0 = _triangles[m][0];
                Vector3<T> v1 = _triangles[m][1];
                Vector3<T> v2 = _triangles[m][2];

                geometry::Triangle<T> CurTri = geometry::Triangle<T>(v0, v1, v2);

                Vector3<T> pmin = v2.cwiseMin(v0.cwiseMin(v1));
                Vector3<T> pmax = v2.cwiseMax(v0.cwiseMax(v1));

                int nxmin = static_cast<int> ((pmin[0] - _pmin[0] + 0.5*_dx)/_dx);
                int nymin = static_cast<int> ((pmin[1] - _pmin[1] + 0.5*_dx)/_dx);

                int nxmax = static_cast<int> ((pmax[0] - _pmin[0]- 0.5*_dx)/_dx);
                int nymax = static_cast<int> ((pmax[1] - _pmin[1]- 0.5*_dx)/_dx);

                for(int i = nxmin; i < nxmax + 1; ++i)
                    for(int j = nymin; j < nymax + 1; ++j)
                    {
                        Vector3<T> base(_pmin[0] + (i + 0.5) * _dx, _pmin[1] + (j + 0.5) * _dx, _pmin[2]);
                        Vector3<T> dir(0,0,1);
                        const T t = CurTri.IntersectRay(base, dir);

                        if (t == -1) continue;
                        _hits[i][j].push_back(base[2] + t * dir[2]);
                    }
            }

			// Step 2: fill the _voxels array
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    {
                        if(_hits[i][j].size()==1) continue;

                        std::sort(_hits[i][j].begin(),_hits[i][j].end(), [](const T A, const T B) -> bool
                        {if (B - A > 1e-7) return true; return false;});

                        for (int idxIn = 1; idxIn < _hits[i][j].size(); idxIn += 2)
                        {
                            int idxOut = idxIn + 1;

                            T zmin = _hits[i][j][idxIn];
                            T zmax = _hits[i][j][idxOut];

                            int nzmin = static_cast<int> ((zmin - _pmin[2] + 0.5*_dx)/_dx);
                            int nzmax = static_cast<int> ((zmax - _pmin[2]- 0.5*_dx)/_dx);

                            for(int n = nzmin; n < nzmax; ++n)
                            {
                                _voxels[i][j][n] = true;
                            }
                        }
                    }
		}

		// TODO: HW2
		// part 3.1.
		// Fill the _voxels array with the correct flag.
		void AdvancedVoxelizationWithApproximation() {
			const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2]; // number of voxels in each dimension
			for (int i = 0; i < nx; ++i)
				for (int j = 0; j < ny; ++j)
					for (int k = 0; k < nz; ++k)
						_voxels[i][j][k] = false;

            std::vector<std::vector<std::vector<int>>> stats = std::vector<std::vector<std::vector<int>>>(nx,
                                                                  std::vector<std::vector<int>>(ny,
                                                                 std::vector<int>(nz, 0)));

            Vector3<T> _pmax = _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx;
            int k = 10;
            for(int i = 0; i < k; ++i)
            {
                // Step 2: randomly generate a direction and compute inside/outside
                Vector3<T> ProjVec = VectorGenerator(i);
                Vector3<T> newmin;
                Vector3<T> newmax;

                ComputerNewBasePlane(ProjVec, _pmin, _pmax, newmin, newmax);
                //std::cout<<ProjVec<<newmin<<newmax;

                const T dx = abs(ProjVec[0])/abs(ProjVec[2])*_dx;
                const T dy = abs(ProjVec[1])/abs(ProjVec[2])*_dx;

                int cnx;
                int cny;
                if(dx ==0)
                {
                    cnx = nx;
                    cny = ny;}
                else{
                cnx = static_cast<int>(newmax[0] - newmin[0]/dx);
                cny = static_cast<int>(newmax[1] - newmin[1]/dy);
                }

                std::vector<std::vector<std::vector<T>>> _hits =
                        std::vector<std::vector<std::vector<T>>>(cnx,
                         std::vector<std::vector<T>>(cny,
                        std::vector<T>(1, -1)));

                //Use ray shooting algorithm from section 2.2
                for(int m = 0; m < _triangles.size(); ++m)
                {
                    Vector3<T> v0 = _triangles[m][0];
                    Vector3<T> v1 = _triangles[m][1];
                    Vector3<T> v2 = _triangles[m][2];

                    Vector3<T> norm(0,0,1);

                    Vector3<T> projv0 = v0 + ProjVec*((newmin - v0).dot(norm))/(norm.dot(ProjVec));
                    Vector3<T> projv1 = v1 + ProjVec*((newmin - v1).dot(norm))/(norm.dot(ProjVec));
                    Vector3<T> projv2 = v2 + ProjVec*((newmin - v2).dot(norm))/(norm.dot(ProjVec));

                    geometry::Triangle<T> CurTri = geometry::Triangle<T>(v0, v1, v2);

                    Vector3<T> trpmin = projv2.cwiseMin(projv0.cwiseMin(projv1));
                    Vector3<T> trpmax = projv2.cwiseMax(projv0.cwiseMax(projv1));

                    int nxmin = static_cast<int> ((trpmin[0] - newmin[0] + 0.5*dx)/dx);
                    int nymin = static_cast<int> ((trpmin[1] - newmin[1] + 0.5*dy)/dy);

                    int nxmax = static_cast<int> ((trpmax[0] - newmin[0]- 0.5*dx)/dx);
                    int nymax = static_cast<int> ((trpmax[1] - newmin[1]- 0.5*dy)/dy);

                    for(int i = nxmin; i < nxmax + 1; ++i)
                        for(int j = nymin; j < nymax + 1; ++j)
                        {
                            Vector3<T> base(newmin[0] + (i + 0.5) * dx, newmin[1] + (j + 0.5) * dy, newmin[2]);
                            const T t = CurTri.IntersectRay(base, ProjVec);

                            if (t == -1) continue;
                            _hits[i][j].push_back(t);
                        }
                }

                // Step 3: fill the _voxels array
                for (int i = 0; i < cnx; ++i)
                    for (int j = 0; j < cny; ++j)
                    {
                        if(_hits[i][j].size()<3) continue;

                        std::sort(_hits[i][j].begin(),_hits[i][j].end(), [](const T A, const T B) -> bool
                        {if (B - A > 1e-7) return true; return false;});

                        for (int idxIn = 1; idxIn < _hits[i][j].size() - 1; idxIn += 2)
                        {
                            int idxOut = idxIn + 1;

                            T tmin = _hits[i][j][idxIn];
                            T tmax = _hits[i][j][idxOut];

                            Vector3<T> ptMin = Vector3<T>(newmin[0] + (i + 0.5) * dx, newmin[1] + (j + 0.5) * dy, newmin[2])
                                    + tmin * ProjVec;
                            Vector3<T> ptMax = Vector3<T>(newmin[0] + (i + 0.5) * dx, newmin[1] + (j + 0.5) * dy, newmin[2])
                                    + tmax * ProjVec;

                            int oxmin = static_cast<int> ((ptMin[0] - _pmin[0])/_dx);
                            int oymin = static_cast<int> ((ptMin[1] - _pmin[1])/_dx);
                            int ozmin = static_cast<int> ((ptMin[2] - _pmin[2])/_dx);

                            int oxmax = static_cast<int> ((ptMax[0] - _pmin[0])/_dx);
                            int oymax = static_cast<int> ((ptMax[1] - _pmin[1])/_dx);
                            int ozmax = static_cast<int> ((ptMax[2] - _pmin[2])/_dx);

                            int dist = oxmax - oxmin;

                            for(int n = 0; n < dist; ++n)
                            {
                                int idx = static_cast<int>(oxmin+n*ProjVec[0]);
                                int idy = static_cast<int>(oymin+n*ProjVec[1]);
                                int idz = static_cast<int>(ozmin+n*ProjVec[2]);
                                stats[idx][idy][idz]++;
                            }
                        }
                    }


            }

            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        if(stats[i][j][k] >= 2)
                        _voxels[i][j][k] = true;
		}

		//Helper Function to compute new base plane
		void ComputerNewBasePlane(Vector3<T> ProjVec, Vector3<T> pmin, Vector3<T> pmax,
                            Vector3<T>& newmin, Vector3<T>& newmax){
            Vector3<T> source;
            Vector3<T> target;
            Vector3<T> n(0,0,1);

            if(ProjVec[2] >0) {
                source = pmax;
                target = pmin;
            }
            else{
                source = pmin;
                target = pmax;
            }

            newmin = Vector3<T>(pmin[0],pmin[1],target[2]);
            newmax = Vector3<T>(pmax[0],pmax[1],target[2]);

            Vector3<T> projMPt;
            std::vector<Vector3<T>> planelist;
            planelist.push_back(source);
            planelist.push_back(target);


		    for(int i = 0; i < 2; ++i)
		        for(int j = 0; j < 2; ++j)
		        {
                    Vector3<T> sPt = Vector3<T>(planelist[i][0],planelist[j][1],source[2]);
                    projMPt = sPt + ProjVec*((target - sPt).dot(n))/(n.dot(ProjVec));

                    newmin = newmin.cwiseMin(projMPt);
                    newmax = newmax.cwiseMax(projMPt);
		        }
		}

		Vector3<T> VectorGenerator(int count)
		{
		    if(count == 8) return Vector3<T>(0,0,1);
            if(count == 9) return Vector3<T>(0,0,-1);
		    Vector3<T> ProjVec(1,1,1);
		    int i = 0;
		    while(count >= 1)
		    {
		        ProjVec[2-i] = pow(-1,(count %2));
		        count = count /2;
		        ++i;
		    }
		    return ProjVec;
		}

		void WriteVoxelToMesh(const std::string& stl_file_name) const {
			const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
			std::vector<std::vector<Vector3<int>>> faces;
			std::vector<Vector3<int>> corners({
				Vector3<int>(0, 0, 0),
				Vector3<int>(0, 0, 1),
				Vector3<int>(0, 1, 0),
				Vector3<int>(0, 1, 1),
				Vector3<int>(1, 0, 0),
				Vector3<int>(1, 0, 1),
				Vector3<int>(1, 1, 0),
				Vector3<int>(1, 1, 1)
				});
			for (int i = 0; i < nx; ++i)
				for (int j = 0; j < ny; ++j)
					for (int k = 0; k < nz; ++k) {
						if (!_voxels[i][j][k]) continue;
						// Check -x direction.
						Vector3<int> cmin(i, j, k);
						if (i == 0 || !_voxels[i - 1][j][k]) {
							faces.push_back({ cmin + corners[0], cmin + corners[1], cmin + corners[3] });
							faces.push_back({ cmin + corners[0], cmin + corners[3], cmin + corners[2] });
						}
						if (i == nx - 1 || !_voxels[i + 1][j][k]) {
							faces.push_back({ cmin + corners[4], cmin + corners[6], cmin + corners[7] });
							faces.push_back({ cmin + corners[4], cmin + corners[7], cmin + corners[5] });
						}
						if (j == 0 || !_voxels[i][j - 1][k]) {
							faces.push_back({ cmin + corners[0], cmin + corners[4], cmin + corners[5] });
							faces.push_back({ cmin + corners[0], cmin + corners[5], cmin + corners[1] });
						}
						if (j == ny - 1 || !_voxels[i][j + 1][k]) {
							faces.push_back({ cmin + corners[2], cmin + corners[3], cmin + corners[7] });
							faces.push_back({ cmin + corners[2], cmin + corners[7], cmin + corners[6] });
						}
						if (k == 0 || !_voxels[i][j][k - 1]) {
							faces.push_back({ cmin + corners[0], cmin + corners[2], cmin + corners[6] });
							faces.push_back({ cmin + corners[0], cmin + corners[6], cmin + corners[4] });
						}
						if (k == nz - 1 || !_voxels[i][j][k + 1]) {
							faces.push_back({ cmin + corners[5], cmin + corners[7], cmin + corners[3] });
							faces.push_back({ cmin + corners[5], cmin + corners[3], cmin + corners[1] });
						}
					}
			std::ofstream fout(stl_file_name);
			fout << "solid vcg" << std::endl;
			for (const auto& f : faces) {
				std::vector<Vector3<T>> p;
				for (const auto& fi : f) {
					Vector3<T> v = _pmin + fi.cast<T>() * _dx;
					p.push_back(v);
				}
				const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
				fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
				fout << "    outer loop" << std::endl;
				for (const auto& v : p) {
					fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
				}
				fout << "    endloop" << std::endl;
				fout << "  endfacet" << std::endl;
			}
			fout << "endsolid vcg" << std::endl;
		}

		void WriteVoxelToFile(const std::string &voxel_file) const {
			// File format: _pmin.x _pmin.y _pmin.z dx nx ny nz
			// Then a [nx+1][ny+1][nz+1] array of 0s and 1s.
			const int nx = _nvoxel.x(), ny = _nvoxel.y(), nz = _nvoxel.z();


			// Write it to file.
			std::ofstream fout(voxel_file);
			fout << _pmin.x() << " " << _pmin.y() << " " << _pmin.z() << " " << _dx
				<< " " << nx << " " << ny << " " << nz << std::endl;

			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					for (int k = 0; k < nz; ++k) {
						if (!_voxels[i][j][k]) {
							fout << "0";
						}
						else {
							fout << "1";
						}
					}
					fout << "0" << std::endl;
				}
				fout << std::string(nz + 1, '0') << std::endl;
			}
			for (int j = 0; j < ny + 1; ++j) {
				fout << std::string(nz + 1, '0') << std::endl;
			}
		}

		using Vector8i = Eigen::Matrix<int, 8, 1>;

		const materials::HexahedralMesh<T> ConvertToHexMesh(std::vector<int>& force, std::vector<int>& fixed, double& num_voxels) {
			std::vector<Vector3<T>> vertices;
			std::vector<Vector8i> elements;
			std::map<int, int> idx_map;
			idx_map.clear();
			// number of possible vertices, not voxels !!
			const int nx = _nvoxel[0] + 1, ny = _nvoxel[1] + 1, nz = _nvoxel[2] + 1;
			int count = 0;
			num_voxels = 0;

			printf("nx: %d, ny: %d, nz: %d\n", nx - 1, ny - 1, nz - 1);

			std::vector<Vector3<int>> corners({
				Vector3<int>(0, 0, 0),
				Vector3<int>(0, 0, 1),
				Vector3<int>(0, 1, 0),
				Vector3<int>(0, 1, 1),
				Vector3<int>(1, 0, 0),
				Vector3<int>(1, 0, 1),
				Vector3<int>(1, 1, 0),
				Vector3<int>(1, 1, 1)
				});

			// loop through number of voxels 
			for (int i = 0; i < nx - 1; ++i)
				for (int j = 0; j < ny - 1; ++j)
					for (int k = 0; k < nz - 1; ++k) {
						if (!_voxels[i][j][k]) continue;
						num_voxels++;

						Vector8i element;
						// push eight vertices to the vertices and elements
						for (int l = 0; l < 8; ++l) {
							int idx = (i + corners[l][0]) * ny * nz + (j + corners[l][1])
								* nz + (k + corners[l][2]);

							if (idx_map.find(idx) == idx_map.end()) {
								idx_map[idx] = count;
								count += 1;
								element(l) = idx_map[idx];
								Vector3<int> vIndex(i + corners[l][0], j + corners[l][1], k + corners[l][2]);
								vertices.push_back(_pmin + vIndex.cast<T>() * _dx);
							}
							else {
								element(l) = idx_map[idx];
							}
						}
						elements.push_back(element);
					}

			// create vertice and element matrix
			materials::Matrix3X<T> vMatrix = MatrixX<T>::Zero(3, vertices.size());
			materials::Matrix8Xi<T> eMatrix = MatrixX<int>::Zero(8, elements.size());

			for (int i = 0; i < vertices.size(); ++i) {
				vMatrix.col(i) = vertices[i];
			}

			for (int i = 0; i < elements.size(); ++i) {
				eMatrix.col(i) = elements[i];
			}

			// find top-z vertices
			for (int i = 0; i < nx - 1; ++i)
				for (int j = 0; j < ny - 1; ++j)
					for (int k = nz - 2; k >= 0; --k) {
						if (!_voxels[i][j][k]) continue;
						printf("i: %d, j: %d, k: %d\n", i, j, k);

						// only store the top points
						for (int l = 1; l < 8; l = l + 2) {
							int idx = (i + corners[l][0]) * ny * nz + (j + corners[l][1])
								* nz + (k + corners[l][2]);
							force.push_back(idx_map[idx]);
						}
						break;
					}
			for (auto i : force) {
				std::cout << i << std::endl;
			}
			std::set<int> forceSet(force.begin(), force.end());
			force.assign(forceSet.begin(), forceSet.end());
			printf("force size: %d\n", force.size());

			std::cout << "top-z done" << std::endl;

			// find left-most vertices
			for (int j = 0; j < ny - 1; ++j)
				for (int k = 0; k < nz - 1; ++k)
					for (int i = 0; i < nx - 1; ++i) {
						if (!_voxels[i][j][k]) continue;

						printf("i: %d, j: %d, k: %d\n", i, j, k);

						// only store the top points
						for (int l = 0; l < 4; l++) {
							int idx = (i + corners[l][0]) * ny * nz + (j + corners[l][1])
								* nz + (k + corners[l][2]);
							fixed.push_back(idx_map[idx]);
						}
						break;
					}
			std::cout << "left-most done" << std::endl;

			// find right-most vertices
			for (int j = 0; j < ny - 1; ++j)
				for (int k = 0; k < nz - 1; ++k)
					for (int i = nx - 2; i >= 0; --i) {
						if (!_voxels[i][j][k]) continue;

						// only store the top points
						for (int l = 4; l < 8; l++) {
							int idx = (i + corners[l][0]) * ny * nz + (j + corners[l][1])
								* nz + (k + corners[l][2]);
							fixed.push_back(idx_map[idx]);
						}
						break;
					}
			std::cout << "right-most done" << std::endl;

			std::set<int> fixedSet(fixed.begin(), fixed.end());
			fixed.assign(fixedSet.begin(), fixedSet.end());
			printf("fixed size: %d\n", fixed.size());

			return materials::HexahedralMesh<T>(vMatrix, eMatrix);
		}

	private:
		std::vector<std::vector<Vector3<T>>> _triangles;
		T _dx;  // The size of each voxel.
		Vector3<T> _pmin;    // The min and max corner of the bounding box.
		Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.
		std::vector<std::vector<std::vector<bool>>> _voxels;   // True <-> voxel is occupied.
	};

}
