#ifndef CGAL_POLYGON_MESH_PROCESSING_LIFTLAYER_MESH_OPERATIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_LIFTLAYER_MESH_OPERATIONS_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <vector>
#include <string>
#include <cmath>

namespace CGAL {
namespace Polygon_mesh_processing {

namespace PMP = CGAL::Polygon_mesh_processing;

/**
 * Structure to hold information about a mesh part
 */
template<typename Kernel>
struct PartInfo {
  Surface_mesh<typename Kernel::Point_3> mesh;
  double avg_z;
  typename Kernel::Point_3 centroid;
  int source; // 1 for critical_mesh, 2 for cut_fill_mesh
  size_t index;
};

/**
 * Compute average Z coordinate of a mesh
 */
template<typename Mesh>
double compute_avg_z(const Mesh& m) {
  auto vpm = get(boost::vertex_point, m);
  double sum_z = 0;
  int count = 0;
  for(auto v : vertices(m)) {
    sum_z += vpm[v].z();
    count++;
  }
  return (count > 0) ? sum_z / count : 0.0;
}

/**
 * Compute centroid of a mesh
 */
template<typename Mesh, typename Kernel>
typename Kernel::Point_3 compute_centroid(const Mesh& m) {
  auto vpm = get(boost::vertex_point, m);
  typename Kernel::Point_3 sum(0, 0, 0);
  int count = 0;
  for(auto v : vertices(m)) {
    auto p = vpm[v];
    sum = typename Kernel::Point_3(sum.x() + p.x(), sum.y() + p.y(), sum.z() + p.z());
    count++;
  }
  return typename Kernel::Point_3(sum.x() / count, sum.y() / count, sum.z() / count);
}

/**
 * Create a cutting surface from the intersection edge between two meshes.
 * Returns true if edge intersection is detected and cutting surface created, false otherwise.
 * This helps split meshes that only intersect at edges.
 * 
 * @param mesh1 First mesh
 * @param mesh2 Second mesh
 * @param cutting_surface Output cutting surface mesh
 * @return true if cutting surface was created, false otherwise
 */
template<typename Mesh, typename Kernel>
bool create_cutting_surface_from_intersection(const Mesh& mesh1, const Mesh& mesh2, Mesh& cutting_surface) {
  // Find the intersection edge(s) using corefine
  Mesh m1_copy = mesh1;
  Mesh m2_copy = mesh2;
  
  auto ecm = m1_copy.template add_property_map<typename boost::graph_traits<Mesh>::edge_descriptor, bool>("e:is_constrained", false).first;
  PMP::corefine(m1_copy, m2_copy, CGAL::parameters::edge_is_constrained_map(ecm), CGAL::parameters::default_values());
  
  // Find constrained edges
  std::vector<typename boost::graph_traits<Mesh>::edge_descriptor> intersection_edges;
  for(auto e : edges(m1_copy)) {
    if(ecm[e]) {
      intersection_edges.push_back(e);
    }
  }
  
  if(intersection_edges.empty()) {
    return false;
  }
  
  std::cout << "Found " << intersection_edges.size() << " edge intersection(s), creating cutting surface\n";
  
  auto vpm = get(boost::vertex_point, m1_copy);
  
  // Use the first intersection edge to define a cutting plane
  auto ie = intersection_edges[0];
  auto v1 = source(ie, m1_copy);
  auto v2 = target(ie, m1_copy);
  auto p1 = vpm[v1];
  auto p2 = vpm[v2];
  
  // Edge direction
  typename Kernel::Vector_3 edge_dir(p1, p2);
  
  // For planar meshes, use Z-axis to find perpendicular direction in XY plane
  typename Kernel::Vector_3 z_axis(0, 0, 1);
  typename Kernel::Vector_3 plane_normal = CGAL::cross_product(edge_dir, z_axis);
  double normal_len = std::sqrt(plane_normal.squared_length());
  
  if(normal_len < 1e-10) {
    // Edge is vertical or near-vertical, use different normal
    plane_normal = typename Kernel::Vector_3(edge_dir.y(), -edge_dir.x(), 0);
    normal_len = std::sqrt(plane_normal.squared_length());
  }
  
  if(normal_len < 1e-10) {
    return false;
  }
  
  // Normalize the plane normal
  plane_normal = plane_normal / normal_len;
  
  // Calculate extent from both meshes to create a large enough cutting surface
  double max_extent = 0.0;
  auto vpm1 = get(boost::vertex_point, mesh1);
  auto vpm2 = get(boost::vertex_point, mesh2);
  
  for(auto v : vertices(mesh1)) {
    auto p = vpm1[v];
    double dist = std::sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
    max_extent = std::max(max_extent, dist);
  }
  for(auto v : vertices(mesh2)) {
    auto p = vpm2[v];
    double dist = std::sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
    max_extent = std::max(max_extent, dist);
  }
  double extension = max_extent * 2.0;
  
  // Create a large rectangular cutting surface perpendicular to the plane normal
  typename Kernel::Point_3 center = CGAL::midpoint(p1, p2);
  
  // Create 4 corners of the cutting rectangle
  typename Kernel::Vector_3 edge_dir_norm = edge_dir / std::sqrt(edge_dir.squared_length());
  typename Kernel::Point_3 c1 = center + edge_dir_norm * extension + plane_normal * extension;
  typename Kernel::Point_3 c2 = center + edge_dir_norm * extension - plane_normal * extension;
  typename Kernel::Point_3 c3 = center - edge_dir_norm * extension - plane_normal * extension;
  typename Kernel::Point_3 c4 = center - edge_dir_norm * extension + plane_normal * extension;
  
  // Create cutting surface mesh
  cutting_surface.clear();
  auto cv1 = cutting_surface.add_vertex(c1);
  auto cv2 = cutting_surface.add_vertex(c2);
  auto cv3 = cutting_surface.add_vertex(c3);
  auto cv4 = cutting_surface.add_vertex(c4);
  
  cutting_surface.add_face(cv1, cv2, cv3);
  cutting_surface.add_face(cv1, cv3, cv4);
  
  std::cout << "Created reusable cutting surface\n";
  
  return true;
}

/**
 * Split two surface meshes and select parts based on Z coordinate (top or bottom).
 * 
 * @param critical_mesh Critical mesh (e.g., existing terrain/surface)
 * @param cut_fill_mesh Cut/fill mesh (e.g., proposed changes)
 * @param select_top If true, select top parts (higher Z); if false, select bottom parts (lower Z)
 * @param output_parts_critical Output vector to store selected parts from critical_mesh
 * @param output_parts_cut_fill Output vector to store selected parts from cut_fill_mesh
 * @param output_prefix Prefix for output files (e.g., "top" or "bottom")
 * 
 * The function splits both meshes by their intersection and selects parts based on Z values.
 * When parts overlap in XY plane and are at same Z level, critical_mesh parts are preferred (face intersection).
 */
template<typename Mesh, typename Kernel>
void split_and_select_surface_meshes(
    const Mesh& critical_mesh, 
    const Mesh& cut_fill_mesh, 
    bool select_top,
    std::vector<Mesh>& output_parts_critical,
    std::vector<Mesh>& output_parts_cut_fill,
    const std::string& output_prefix = "selected")
{
  std::cout << "\n=== Collecting all parts from both meshes ===\n";
  
  // Try to create a cutting surface for edge-only intersections
  Mesh cutting_surface;
  bool has_cutting_surface = create_cutting_surface_from_intersection<Mesh, Kernel>(critical_mesh, cut_fill_mesh, cutting_surface);
  
  // Split critical_mesh
  std::vector<Mesh> parts_critical;
  Mesh critical_split = critical_mesh;
  
  try {
    if (has_cutting_surface) {
      std::cout << "Using cutting surface to split critical_mesh\n";
      PMP::split(critical_split, cutting_surface);
    } else {
      std::cout << "Using standard split for critical_mesh\n";
      Mesh cut_fill_for_split = cut_fill_mesh;
      PMP::split(critical_split, cut_fill_for_split);
    }
    
    typedef typename boost::graph_traits<Mesh>::faces_size_type faces_size_type;
    auto pidmap_critical = critical_split.template add_property_map<typename boost::graph_traits<Mesh>::face_descriptor, faces_size_type>("f:patch_id_critical", 0).first;
    PMP::connected_components(critical_split, pidmap_critical);
    PMP::split_connected_components(critical_split, parts_critical, CGAL::parameters::face_patch_map(pidmap_critical));
    std::cout << "Critical mesh split into " << parts_critical.size() << " parts\n";
  } catch (const std::exception& e) {
    std::cerr << "Error splitting critical_mesh: " << e.what() << std::endl;
    parts_critical.push_back(critical_mesh);
  }
  
  // Split cut_fill_mesh using the same cutting surface
  std::vector<Mesh> parts_cut_fill;
  Mesh cut_fill_split = cut_fill_mesh;
  
  try {
    if (has_cutting_surface) {
      std::cout << "Using cutting surface to split cut_fill_mesh\n";
      PMP::split(cut_fill_split, cutting_surface);
    } else {
      std::cout << "Using standard split for cut_fill_mesh\n";
      Mesh critical_for_split = critical_mesh;
      PMP::split(cut_fill_split, critical_for_split);
    }
    
    typedef typename boost::graph_traits<Mesh>::faces_size_type faces_size_type;
    auto pidmap_cut_fill = cut_fill_split.template add_property_map<typename boost::graph_traits<Mesh>::face_descriptor, faces_size_type>("f:patch_id_cut_fill", 0).first;
    PMP::connected_components(cut_fill_split, pidmap_cut_fill);
    PMP::split_connected_components(cut_fill_split, parts_cut_fill, CGAL::parameters::face_patch_map(pidmap_cut_fill));
    std::cout << "Cut/fill mesh split into " << parts_cut_fill.size() << " parts\n";
  } catch (const std::exception& e) {
    std::cerr << "Error splitting cut_fill_mesh: " << e.what() << std::endl;
    parts_cut_fill.push_back(cut_fill_mesh);
  }
  
  // Store all parts with their metadata
  std::vector<PartInfo<Kernel>> all_parts;
  
  for(size_t i = 0; i < parts_critical.size(); ++i) {
    PartInfo<Kernel> info;
    info.mesh = parts_critical[i];
    info.avg_z = compute_avg_z(parts_critical[i]);
    info.centroid = compute_centroid<Mesh, Kernel>(parts_critical[i]);
    info.source = 1; // 1 for critical_mesh
    info.index = i;
    all_parts.push_back(info);
  }
  
  for(size_t i = 0; i < parts_cut_fill.size(); ++i) {
    PartInfo<Kernel> info;
    info.mesh = parts_cut_fill[i];
    info.avg_z = compute_avg_z(parts_cut_fill[i]);
    info.centroid = compute_centroid<Mesh, Kernel>(parts_cut_fill[i]);
    info.source = 2; // 2 for cut_fill_mesh
    info.index = i;
    all_parts.push_back(info);
  }
  
  std::cout << "\n=== Part information ===\n";
  for(size_t i = 0; i < all_parts.size(); ++i) {
    std::string source_name = (all_parts[i].source == 1) ? "critical" : "cut_fill";
    std::cout << "Part " << i << " (from " << source_name << "_mesh): "
              << "avg_z = " << all_parts[i].avg_z 
              << ", centroid = (" << all_parts[i].centroid.x() << ", " 
              << all_parts[i].centroid.y() << ", " << all_parts[i].centroid.z() << ")\n";
  }
  
  // Selection strategy based on select_top parameter
  const double z_tolerance = 0.05; // Tolerance to consider parts at same Z level
  const double xy_overlap_threshold = 5.0; // Distance threshold for XY plane overlap detection
  std::vector<size_t> selected_indices;
  
  std::cout << "Selecting " << (select_top ? "TOP" : "BOTTOM") << " parts\n";

  //TODO the selection logic can be improved. Given a different size of map,the xy_overlay_threshold might need to be adjusted accordingly.
  //Otherwise, the logic might fail to detect overlapping parts in some cases.
  
  // Simple approach: For each part, check if there's any other part that's 
  // higher/lower (depending on mode) and overlapping in XY plane
  for(size_t i = 0; i < all_parts.size(); ++i) {
    bool is_topmost = true;  // For top mode: assume this part is on top until proven otherwise
    
    for(size_t j = 0; j < all_parts.size(); ++j) {
      if(i == j) continue;
      
      // Check if parts overlap in XY by comparing centroids
      double dx = all_parts[i].centroid.x() - all_parts[j].centroid.x();
      double dy = all_parts[i].centroid.y() - all_parts[j].centroid.y();
      double dist_xy = std::sqrt(dx*dx + dy*dy);
      
      // Debug: Print distance information
      std::cout << "  Comparing Part " << i << " vs Part " << j 
                << ": XY distance = " << dist_xy 
                << " (threshold = " << xy_overlap_threshold << ")";
      
      // Consider parts overlapping if their centroids are close
      bool overlaps = (dist_xy < xy_overlap_threshold);
      
      if(overlaps) {
        std::cout << " -> OVERLAPS\n";
        if(select_top) {
          // For top mode: exclude if another part is higher
          if(all_parts[j].avg_z > all_parts[i].avg_z + z_tolerance) {
            is_topmost = false;
            std::cout << "Part " << i << " (" << (all_parts[i].source == 1 ? "critical" : "cut_fill") << ", z=" << all_parts[i].avg_z 
                      << ") excluded: covered by part " << j << " (" << (all_parts[j].source == 1 ? "critical" : "cut_fill") 
                      << ", z=" << all_parts[j].avg_z << ")\n";
            break;
          }
          // If at same Z level (within tolerance), prefer critical_mesh
          else if(std::abs(all_parts[j].avg_z - all_parts[i].avg_z) <= z_tolerance) {
            if(all_parts[j].source == 1 && all_parts[i].source == 2) {
              is_topmost = false;
              std::cout << "Part " << i << " (cut_fill, z=" << all_parts[i].avg_z 
                        << ") excluded: same level as critical part " << j << " (face intersection)\n";
              break;
            }
          }
        } else {
          // For bottom mode: exclude if another part is lower
          if(all_parts[j].avg_z < all_parts[i].avg_z - z_tolerance) {
            is_topmost = false;
            std::cout << "Part " << i << " (" << (all_parts[i].source == 1 ? "critical" : "cut_fill") << ", z=" << all_parts[i].avg_z 
                      << ") excluded: covered by part " << j << " (" << (all_parts[j].source == 1 ? "critical" : "cut_fill") 
                      << ", z=" << all_parts[j].avg_z << ")\n";
            break;
          }
          // If at same Z level (within tolerance), prefer cut_fill_mesh for bottom mode
          else if(std::abs(all_parts[j].avg_z - all_parts[i].avg_z) <= z_tolerance) {
            if(all_parts[j].source == 2 && all_parts[i].source == 1) {
              is_topmost = false;
              std::cout << "Part " << i << " (critical, z=" << all_parts[i].avg_z 
                        << ") excluded: same level as cut_fill part " << j << " (face intersection)\n";
              break;
            }
          }
        }
      } else {
        std::cout << " -> no overlap\n";
      }
    }
    
    if(is_topmost) {
      selected_indices.push_back(i);
    }
  }
  
  std::cout << "\n=== Selected " << (select_top ? "TOP" : "BOTTOM") << " parts ===\n";
  for(size_t idx : selected_indices) {
    std::string source_name = (all_parts[idx].source == 1) ? "critical" : "cut_fill";
    std::cout << "Selected part " << idx << " from " << source_name << "_mesh"
              << " (avg_z = " << all_parts[idx].avg_z << ")\n";
  }
  
  // Separate selected parts by source and store in output vectors
  output_parts_critical.clear();
  output_parts_cut_fill.clear();
  
  int critical_count = 0;
  int cut_fill_count = 0;
  
  for(size_t idx : selected_indices) {
    if(all_parts[idx].source == 1) {
      // Critical mesh part
      std::string filename = output_prefix + "_part_" + std::to_string(critical_count) + "_from_critical.off";
      CGAL::IO::write_polygon_mesh(filename, all_parts[idx].mesh, CGAL::parameters::stream_precision(17));
      std::cout << output_prefix << " part " << critical_count << " from critical_mesh written to " << filename << "\n";
      output_parts_critical.push_back(all_parts[idx].mesh);
      critical_count++;
    } else {
      // Cut/fill mesh part
      std::string filename = output_prefix + "_part_" + std::to_string(cut_fill_count) + "_from_cut_fill.off";
      CGAL::IO::write_polygon_mesh(filename, all_parts[idx].mesh, CGAL::parameters::stream_precision(17));
      std::cout << output_prefix << " part " << cut_fill_count << " from cut_fill_mesh written to " << filename << "\n";
      output_parts_cut_fill.push_back(all_parts[idx].mesh);
      cut_fill_count++;
    }
  }
  
  std::cout << "Total selected: " << critical_count << " critical parts, " 
            << cut_fill_count << " cut/fill parts\n";
  
  // Write all parts for debugging
  std::cout << "\n=== Writing all individual parts for inspection ===\n";
  for(size_t i = 0; i < all_parts.size(); ++i) {
    std::string source_name = (all_parts[i].source == 1) ? "critical" : "cut_fill";
    std::string filename = "all_part_" + std::to_string(i) + "_" + source_name + ".off";
    CGAL::IO::write_polygon_mesh(filename, all_parts[i].mesh, CGAL::parameters::stream_precision(17));
  }
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_LIFTLAYER_MESH_OPERATIONS_H
