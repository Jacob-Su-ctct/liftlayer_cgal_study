#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
//#include <CGAL/Polygon_mesh_processing/barycentric_coordinates.h> // TODO
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>


#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                                          Point_3;
typedef K::Point_2                                          Point_2;
typedef K::Segment_3                                        Segment_3;
typedef K::Ray_3                                            Ray_3;
typedef K::Vector_3                                         Vector_3;
typedef CGAL::Surface_mesh<Point_3>                         Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh>      Primitive;
typedef CGAL::AABB_traits_3<K, Primitive>                   AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>                        Tree;
typedef K::Aff_transformation_3                             Transform;

typedef SMS::Edge_length_cost<Mesh>   Cost;
typedef SMS::Midpoint_placement<Mesh> Placement;

bool extrude_and_clip(Mesh& toClip, Mesh& clipper, double verticalExtrusion, bool reverseFaces)
{
  if (clipper.number_of_faces() == 0)
  {
    return false;
  }

  bool extrudeDown = verticalExtrusion < 0;

  Mesh extrudedClipper;
  PMP::extrude_mesh(clipper, extrudedClipper, Vector_3(0, 0, std::abs(verticalExtrusion)));
  if (reverseFaces) PMP::reverse_face_orientations(extrudedClipper);
  if (extrudeDown)
  {
    Transform translateDown(CGAL::Translation(), Vector_3(0, 0, verticalExtrusion));
    PMP::transform(translateDown, extrudedClipper);
  }
  PMP::clip(toClip, extrudedClipper, CGAL::parameters::do_not_modify(true));

  // Clean up clipped surface
  Mesh temp;
  CGAL::copy_face_graph(toClip, temp);
  toClip.clear();
  CGAL::copy_face_graph(temp, toClip);

  PMP::remove_degenerate_faces(toClip);

  return true;
}

bool is_surface_equal(Mesh& a, Mesh& b)
{
  return a.number_of_faces() == b.number_of_faces() && a.number_of_edges() == b.number_of_edges() && a.number_of_halfedges() == b.number_of_halfedges();
}

void combine_surfaces(Mesh& a, Mesh& b, Mesh& c, Mesh& output)
{
  CGAL::copy_face_graph(a, output);
  CGAL::copy_face_graph(b, output);
  CGAL::copy_face_graph(c, output);

  //PMP::repair_polygon_soup(output, CGAL::parameters::vertex_point_map(output.points()).geom_traits(K()));
  PMP::remove_degenerate_faces(output);
  PMP::remove_isolated_vertices(output);
}

void generate_lift_layers(Mesh& critical_in, Mesh& cut_in, Mesh& fill_in)
{
  Mesh composite;
  int currentLayer = 0;
  double cutOffset = 0.2;
  double fillOffset = 0.35;
  bool has_fill_component = true;
  bool has_cut_component = true;
  double totalTime = 0;

  while (has_fill_component && has_cut_component && currentLayer < 100)
  {
    auto start = std::chrono::high_resolution_clock::now();

    composite.clear();
    Mesh critical = critical_in;
    Mesh cut = cut_in;
    Transform translateDown(CGAL::Translation(), Vector_3(0, 0, -cutOffset * (double)currentLayer));
    PMP::transform(translateDown, cut);
    Mesh fill = fill_in;
    Transform translateUp(CGAL::Translation(), Vector_3(0, 0, fillOffset * (double)currentLayer));
    PMP::transform(translateUp, fill);

    // Step 1: Clip the cut mesh to be only the segments of the original cut surface that make up the final composite layer
    extrude_and_clip(cut, critical, 5, true);

    // Step 2: Clip the fill mesh to only be the segments of the original fill surface that make up the final composite layer
    extrude_and_clip(fill, critical, -5, true);

    // Step 3: Clip the critical surface to exclude all parts that are above the final fill segment(s)
    has_fill_component = extrude_and_clip(critical, fill, 5, false);

    // Step 4: Clip the critical surface to exclude all parts that are below the final cut segment(s)
    has_cut_component = extrude_and_clip(critical, cut, -5, false);

    // Step 5: Union all segments of the final composite surface
    combine_surfaces(critical, cut, fill, composite);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    totalTime += duration.count() / 1000.0;

    std::string name = ("layer-" + std::to_string(currentLayer) + ".off");
    CGAL::IO::write_polygon_mesh(name, composite, CGAL::parameters::stream_precision(17));

    currentLayer++;
  }

  std::cout << "Processing time: " << totalTime << std::endl;
}

int main(int argc, char* argv[])
{
 // Parse command line arguments
  std::string criticalSurfaceFile;
  std::string cutSurfaceFile;
  std::string fillSurfaceFile;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <critical.off> --cut <cut.off> --fill <fill.off>" << std::endl;
    return 1;
  }

  criticalSurfaceFile = argv[1];
  for (int i = 2; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--cut" && i + 1 < argc) {
      cutSurfaceFile = argv[i + 1];
      ++i;
    } else if (arg == "--fill" && i + 1 < argc) {
      fillSurfaceFile = argv[i + 1];
      ++i;
    }
  }

  if (criticalSurfaceFile.empty() || cutSurfaceFile.empty() || fillSurfaceFile.empty()) {
    std::cerr << "Missing required arguments. Usage: " << argv[0] << " <critical.off> --cut <cut.off> --fill <fill.off>" << std::endl;
    return 1;
  }

  Mesh critical, cut, fill;
  if(!PMP::IO::read_polygon_mesh(criticalSurfaceFile, critical) || !PMP::IO::read_polygon_mesh(cutSurfaceFile, cut) || !PMP::IO::read_polygon_mesh(fillSurfaceFile, fill))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  generate_lift_layers(critical, cut, fill);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Total time: " << duration.count() / 1000 << std::endl;

  return 0;
}