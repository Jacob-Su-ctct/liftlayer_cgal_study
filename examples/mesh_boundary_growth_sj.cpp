#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Origin.h>
#include <boost/property_map/property_map.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator vertex_iterator;

int main(int argc, char* argv[])
{
  std::string mesh_file = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/plane.off");
  double offset_distance = (argc > 2) ? std::stod(argv[2]) : 0.1;

  Mesh mesh;
  std::ifstream input(mesh_file);

  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Cannot open " << mesh_file << std::endl;
    return 1;
  }

  // Compute vertex normals
  std::vector<Vector> normals(mesh.number_of_vertices());
  auto vnormals = boost::make_iterator_property_map(normals.begin(), get(boost::vertex_index, mesh));
  CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormals);

  // Compute mesh centroid to determine outward direction
  Vector sum = CGAL::NULL_VECTOR;
  std::size_t vertex_count = 0;
  for (vertex_descriptor v : vertices(mesh)) {
    sum = sum + (mesh.point(v) - CGAL::ORIGIN);
    ++vertex_count;
  }
  Point centroid = CGAL::ORIGIN;
  if (vertex_count > 0) {
    centroid = CGAL::ORIGIN + (sum / static_cast<double>(vertex_count));
  }

  // Move boundary vertices
  for (vertex_descriptor v : vertices(mesh)) {
    bool is_boundary = false;
    Vector boundary_tangent = CGAL::NULL_VECTOR;
    for (auto h : halfedges_around_target(v, mesh)) {
      if (is_border(h, mesh)) {
        is_boundary = true;
        Point src = mesh.point(source(h, mesh));
        Point tgt = mesh.point(target(h, mesh));
        boundary_tangent = boundary_tangent + (src - tgt);
      }
    }
    if (is_boundary) {
      Point current_pos = mesh.point(v);
      Vector normal = normals[v];

      if (normal == CGAL::NULL_VECTOR)
        continue;

      Vector radial = current_pos - centroid;
      Vector offset_dir = radial - (radial * normal) * normal; // project onto tangent plane

      double offset_length = std::sqrt(offset_dir.squared_length());
      if (offset_length == 0.0) {
        if (boundary_tangent != CGAL::NULL_VECTOR) {
          offset_dir = CGAL::cross_product(normal, boundary_tangent);
          offset_length = std::sqrt(offset_dir.squared_length());
        }
      }

      if (offset_length == 0.0)
        continue; // no meaningful direction

      offset_dir = offset_dir / offset_length;

      if (offset_dir * radial < 0.0)
        offset_dir = -offset_dir;

      Point new_pos = current_pos + offset_distance * offset_dir;
      mesh.point(v) = new_pos;
    }
  }

  // Save the offset mesh
  std::ofstream output("offset_mesh.off");
  output << mesh;
  output.close();

  std::cout << "Offset mesh saved to offset_mesh.off" << std::endl;

  return 0;
}