#pragma once
#include <Eigen/Dense>
#include <vector>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    template <typename T>
    using Mat3 = Eigen::Matrix<T, 3, 3>;

    // the Plane rs represented by (x - _p) /dot _normal = 0
    template <typename T>
    class Plane {
    public:
        Plane(Vector3<T> p, Vector3<T> normal) {
            _p = p;
            _normal = normal;
            _normal.normalize();
        }

        Vector3<T>& p() { return _p; }
        Vector3<T>& normal() { return _normal; }
        
        // return if the point is on Plane
        // also fill parameter dist as the signed distance from point to Plane
        bool onPlane(Vector3<T> point, T& dist) {
            dist = (point - _p).dot(_normal);
            if (std::fabs(dist) < 1e-6) {
                return true;
            } else {
                return false;
            }
        }

    private:
        Vector3<T> _p;
        Vector3<T> _normal;
    };

    template <typename T>
    class Triangle {
    public:
        Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;
            _n = (_vertices[2] - _vertices[0]).cross(_vertices[1] - _vertices[0]);
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }

        
        // TODO: HW1
        // part 2.1
        // Implement the function to do intersection between triangle and Plane p
        // Input: Plane p
        // Output: return pairs for intersections with three edges
        // Hint:
        //      - enumerate three edges of the triangle and do intersection individually
        //      - consider the case that no intersection
        //      - consider how to avoid repeated intersection points in returned list
        std::vector<std::pair<int, Vector3<T>>> IntersectPlane(Plane<T> p) {

            // The first entry of the pair marks the intersecting situation
            // -1 - no intersection
            // 0, 1, 2 intersect at vertex 0, 1, 2
            // 3, 4, 5 marks intersection at edge (0,1), (0,2), (1,2)

            std::vector<std::pair<int,Vector3<T>>> intersections;
//            //intersections.clear();
//
//            Vector3<T> n = p.normal();
//            Vector3<T> pt = p.p();
//
//            Vector3<T>& v0 = _vertices[0];
//            Vector3<T>& v1 = _vertices[1];
//            Vector3<T>& v2 = _vertices[2];
//
//            T d0 = n.dot(v0 - pt);
//            T d1 = n.dot(v1 - pt);
//            T d2 = n.dot(v2 - pt);
//
//            Vector3<T> s;
//            Vector3<T> e;
//
//            bool nInt = false;
//
//            if(abs(d0) < 1e-6)
//            {
//                if(abs(d1) < 1e-6)
//                {
//                    if(abs(d2) > 1e-6)
//                    {
//                        intersections.push_back(std::pair<int, Vector3<T>>(0, v0));
//                        intersections.push_back(std::pair<int, Vector3<T>>(1, v1));
//                    }
//                    else nInt = true;
//                }
//                else
//                {
//                    if(abs(d2) < 1e-6)
//                    {
//                        intersections.push_back(std::pair<int, Vector3<T>>(0, v0));
//                        intersections.push_back(std::pair<int, Vector3<T>>(2, v2));
//                    }
//                    else if(d1 * d2 < 0)
//                    {
//                        intersections.push_back(std::pair<int, Vector3<T>>(0, v0));
//                        Vector3<T> v12 = PlaneSegmentIntersection(p, v1, v2);
//                        intersections.push_back(std::pair<int, Vector3<T>>(4, v12));
//                    }
//                    else nInt = true;
//                }
//            }
//            else {
//                if (abs(d1) < 1e-6) {
//                    if (abs(d2) < 1e-6) {
//                        intersections.push_back(std::pair<int, Vector3<T>>(1, v1));
//                        intersections.push_back(std::pair<int, Vector3<T>>(2, v2));
//                    } else if (d0 * d2 < 0) {
//                        intersections.push_back(std::pair<int, Vector3<T>>(1, v1));
//                        // Calculate Intersection of v0v2 and Plane
//                        Vector3<T> v20 = PlaneSegmentIntersection(p, v2, v0);
//                        intersections.push_back(std::pair<int, Vector3<T>>(5, v20));
//                    } else nInt = true;
//                } else {
//                    if (abs(d2) < 1e-6) {
//                        if (d0 * d1 < 0) {
//                            intersections.push_back(std::pair<int, Vector3<T>>(2, v2));
//                            // Calculate Intersection of v0v1
//                            Vector3<T> v01 = PlaneSegmentIntersection(p, v0, v1);
//                            intersections.push_back(std::pair<int, Vector3<T>>(3, v01));
//                        } else nInt = true;
//                    } else {
//                        if (d0 * d2 < - 1e-6 && d1 * d2 < - 1e-6) {
//                            Vector3<T> v20 = PlaneSegmentIntersection(p, v2, v0);
//                            intersections.push_back(std::pair<int, Vector3<T>>(5, v20));
//                            Vector3<T> v12 = PlaneSegmentIntersection(p, v1, v2);
//                            intersections.push_back(std::pair<int, Vector3<T>>(4, v12));
//                        } else if (d0 * d1 < - 1e-6 && d1 * d2 < - 1e-6) {
//                            Vector3<T> v01 = PlaneSegmentIntersection(p, v0, v1);
//                            intersections.push_back(std::pair<int, Vector3<T>>(3, v01));
//                            Vector3<T> v12 = PlaneSegmentIntersection(p, v1, v2);
//                            intersections.push_back(std::pair<int, Vector3<T>>(4, v12));
//                        } else if (d0 * d1 < - 1e-6 && d0 * d2 < - 1e-6) {
//                            Vector3<T> v01 = PlaneSegmentIntersection(p, v0, v1);
//                            intersections.push_back(std::pair<int, Vector3<T>>(3, v01));
//                            Vector3<T> v20 = PlaneSegmentIntersection(p, v0, v2);
//                            intersections.push_back(std::pair<int, Vector3<T>>(5, v20));
//                        } else nInt = true;
//                    }
//                }
//            }
//
//            if(nInt)
//            {
//                intersections.push_back(std::pair<int, Vector3<T>>(-1, s));
//                intersections.push_back(std::pair<int, Vector3<T>>(-1, e));
//            }

            return intersections;
        }

//        Vector3<T> PlaneSegmentIntersection(Plane<T> p, Vector3<T> s, Vector3<T> t)
//        {
//            T lambda = ((p.p() - s).dot(p.normal()))/(t - s).dot(p.normal());
//            Vector3<T> intersection = s + lambda * (t - s);
//            return intersection;
//        }

        // TODO: HW2
        // part 1.1
        // Implement the function to do intersection between triangle and a ray
        // Input: a ray, the ray is represented by an origin position and a direction vector
        // Output: return a real number t, the intersection is origin + dir * t, t = -1 means no intersection
        const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const {

            const T flag = static_cast<T>(-1.0);

            T t;

            const Vector3<T> v0 = _vertices[0];
            const Vector3<T> v1 = _vertices[1];
            const Vector3<T> v2 = _vertices[2];

            if(dir.dot(_n) == 0) return flag;
            else{
                t = (v0 - origin).dot(_n)/(dir.dot(_n));
            }
            Vector3<T> pt = origin + t * dir;

            //http://geomalgorithms.com/a06-_intersect-2.html
            Vector3<T> w = pt - v0;
            Vector3<T> u = v1 - v0;
            Vector3<T> v = v2 - v0;
            const T uv = u.dot(v);
            const T uu = u.dot(u);
            const T vv = v.dot(v);
            const T wu = w.dot(u);
            const T wv = w.dot(v);
            const T s = (uv * wv-vv * wu)/(uv * uv - uu * vv);
            const T r = (uv * wu-uu * wv)/(uv * uv-uu * vv);

            if (s>= 0 &&  r>= 0 && s + r <= 1) return t;
            return flag;

        }

    private:
        Vector3<T> _vertices[3];
        Vector3<T> _n;
    };
}