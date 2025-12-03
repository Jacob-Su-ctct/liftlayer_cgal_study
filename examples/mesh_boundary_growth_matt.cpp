#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <vector>
#include <iostream>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

void extrude_along_normals(Mesh& input_mesh, double offset_distance, Mesh& output_mesh) {
  // 1. Compute Vertex Normals and Offset Positions
  // We need 'angle weighted' normals for better corner behavior
  auto vnormals = input_mesh.add_property_map<vertex_descriptor, K::Vector_3>("v:normals", K::Vector_3(0, 0, 0)).first;
  PMP::compute_vertex_normals(input_mesh, vnormals);

  std::vector<K::Point_3> new_points;
  std::vector<std::vector<std::size_t>> new_polygons;

  // We build the new geometry in a 'soup' container first
  for (vertex_descriptor v : input_mesh.vertices()) {
    K::Point_3 p = input_mesh.point(v);
    K::Vector_3 n = vnormals[v];

    // --- ANGLE CORRECTION (Robust Check) ---
    double scale_factor = 1.0;

    auto h = input_mesh.halfedge(v);

    // 1. Check if vertex is not isolated
    if (h != Mesh::null_halfedge()) {

      // 2. Check if the default halfedge points to a Null Face (Boundary)
      if (input_mesh.face(h) == Mesh::null_face()) {
        // If so, flip to the opposite halfedge, which usually points 'inwards' to a valid face
        h = input_mesh.opposite(h);
      }

      // 3. Only compute normal if we successfully found a valid face
      if (input_mesh.face(h) != Mesh::null_face()) {
        K::Vector_3 fn = PMP::compute_face_normal(input_mesh.face(h), input_mesh);

        double cos_theta = n * fn; // Dot product

        // Protect against division by zero or negative flip
        if (cos_theta > 1e-6) {
          scale_factor = 1.0 / cos_theta;
        }
      }
    }

    // Clamp scaling to avoid exploding spikes at sharp corners
    if (scale_factor > 3.0) scale_factor = 3.0;

    K::Point_3 p_new = p + n * (offset_distance * scale_factor);
    new_points.push_back(p_new);
  }

  // 2. Build the Connectivity (Same as input)
  for (face_descriptor f : input_mesh.faces()) {
    std::vector<std::size_t> poly;
    for (vertex_descriptor v : input_mesh.vertices_around_face(input_mesh.halfedge(f))) {
      poly.push_back(size_t(v));
    }
    new_polygons.push_back(poly);
  }

  // 3. Fix Self-Intersections (The CGAL 6.1 Magic)
  // This splits triangles where they overlap
  PMP::autorefine_triangle_soup(new_points, new_polygons);
  PMP::repair_polygon_soup(new_points, new_polygons);
  PMP::orient_polygon_soup(new_points, new_polygons);

  // 4. Convert back to Mesh
  if (PMP::is_polygon_soup_a_polygon_mesh(new_polygons)) {
    PMP::polygon_soup_to_polygon_mesh(new_points, new_polygons, output_mesh);
  }
  else {
    std::cerr << "Error: The offset resulted in non-manifold geometry (e.g., edges shared by >2 faces)." << std::endl;
    // Fallback or Alpha Wrap required here
  }
}

int main(int argc, char* argv[])
{
  auto now = std::chrono::high_resolution_clock::now();
  auto timeElapsed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());

  const std::string criticalSurfaceFile = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/plane.off");
  double offset_distance = (argc > 2) ? std::stod(argv[2]) : 1.5;

  Mesh critical, cut, fill;
  if(!PMP::IO::read_polygon_mesh(criticalSurfaceFile, critical) )
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  Mesh output;
  extrude_along_normals(critical, offset_distance, output);
  std::string name = ("mesh_extruded_along_normal.off" );
  CGAL::IO::write_polygon_mesh(name, output, CGAL::parameters::stream_precision(17));
  //generate_lift_layers(critical, cut, fill);
  //PMP::clip(critical, cut);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Total time: " << duration.count() / 1000 << std::endl;

  return 0;
}