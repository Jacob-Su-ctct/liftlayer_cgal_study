#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/liftlayer_mesh_operations.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

// Helper function to join multiple meshes into a single composite mesh
Mesh join_meshes(const std::vector<Mesh>& meshes) {
  if(meshes.empty()) {
    return Mesh();
  }
  
  Mesh composite = meshes[0];
  int consecutive_failures = 0;
  const int MAX_CONSECUTIVE_FAILURES = 3;
  
  for(size_t i = 1; i < meshes.size(); ++i) {
    // Check if meshes actually intersect before attempting union
    // Union only works well for overlapping/connected meshes
    bool meshes_touch = false;
    try {
      meshes_touch = PMP::do_intersect(composite, meshes[i]);
    } catch(const std::exception& e) {
      std::cerr << "Warning: Intersection test failed for mesh " << i << ": " << e.what() << "\n";
      meshes_touch = false;
    } catch(...) {
      std::cerr << "Warning: Intersection test failed for mesh " << i << " (unknown error)\n";
      meshes_touch = false;
    }
    
    if(!meshes_touch) {
      // Meshes don't intersect - just append without union
      std::cerr << "Meshes don't intersect (composite vs mesh " << i << "), appending directly.\n";
      
      try {
        std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, 
                 typename boost::graph_traits<Mesh>::vertex_descriptor> vmap;
        auto vpm_src = get(boost::vertex_point, meshes[i]);
        
        for(auto v : vertices(meshes[i])) {
          auto new_v = composite.add_vertex(vpm_src[v]);
          vmap[v] = new_v;
        }
        
        for(auto f : faces(meshes[i])) {
          std::vector<typename boost::graph_traits<Mesh>::vertex_descriptor> face_verts;
          for(auto v : vertices_around_face(halfedge(f, meshes[i]), meshes[i])) {
            face_verts.push_back(vmap[v]);
          }
          if(face_verts.size() >= 3) {
            composite.add_face(face_verts);
          }
        }
      } catch(const std::exception& e) {
        std::cerr << "Error during direct append of mesh " << i << ": " << e.what() << "\n";
      }
      continue;  // Skip union attempt
    }
    
    // Meshes intersect - try union
    // But skip union if the composite is getting too complex (reduces crash risk)
    size_t composite_complexity = num_vertices(composite) + num_faces(composite);
    size_t mesh_complexity = num_vertices(meshes[i]) + num_faces(meshes[i]);
    
    if(composite_complexity > 50 || mesh_complexity > 30) {
      // Complexity threshold exceeded - use simple append to avoid CGAL corefinement issues
      std::cerr << "Mesh complexity high (composite: " << composite_complexity 
                << ", mesh " << i << ": " << mesh_complexity << "), using append instead of union.\n";
      
      try {
        std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, 
                 typename boost::graph_traits<Mesh>::vertex_descriptor> vmap;
        auto vpm_src = get(boost::vertex_point, meshes[i]);
        
        for(auto v : vertices(meshes[i])) {
          auto new_v = composite.add_vertex(vpm_src[v]);
          vmap[v] = new_v;
        }
        
        for(auto f : faces(meshes[i])) {
          std::vector<typename boost::graph_traits<Mesh>::vertex_descriptor> face_verts;
          for(auto v : vertices_around_face(halfedge(f, meshes[i]), meshes[i])) {
            face_verts.push_back(vmap[v]);
          }
          if(face_verts.size() >= 3) {
            composite.add_face(face_verts);
          }
        }
      } catch(const std::exception& e) {
        std::cerr << "Error during complexity-based append of mesh " << i << ": " << e.what() << "\n";
      }
      continue;  // Skip union attempt
    }
    
    std::cerr << "Meshes intersect (composite vs mesh " << i << "), attempting union...\n";
    Mesh temp = composite;
    Mesh current_mesh = meshes[i]; // Create a non-const copy for union operation
    composite.clear();
    
    bool union_succeeded = false;
    try {
      union_succeeded = PMP::corefine_and_compute_union(temp, current_mesh, composite);
      if(union_succeeded) {
        consecutive_failures = 0;  // Reset counter on success
        std::cerr << "Union succeeded for mesh " << i << "\n";
      }
    } catch(const std::exception& e) {
      std::cerr << "Error during union: " << e.what() << "\n";
      union_succeeded = false;
    } catch(...) {
      std::cerr << "Unknown error during union\n";
      union_succeeded = false;
    }
    
    if(!union_succeeded) {
      // If union fails, just append the mesh by copying vertices and faces from original mesh
      std::cerr << "Warning: Union operation failed for mesh " << i << ", appending instead.\n";
      composite = temp; // Start with previous composite
      consecutive_failures++;
      
      // If too many failures, the composite mesh might be getting corrupted
      // In that case, just keep what we have and skip remaining meshes
      if(consecutive_failures >= MAX_CONSECUTIVE_FAILURES) {
        std::cerr << "Too many consecutive union failures (" << consecutive_failures 
                  << "), keeping composite as-is and skipping remaining " 
                  << (meshes.size() - i) << " meshes.\n";
        break;
      }
      
      // Simple append by adding all vertices and faces from the ORIGINAL mesh (meshes[i])
      // Don't use current_mesh as it may have been corrupted by failed union
      std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, 
               typename boost::graph_traits<Mesh>::vertex_descriptor> vmap;
      auto vpm_src = get(boost::vertex_point, meshes[i]);
      auto vpm_tgt = get(boost::vertex_point, composite);
      
      size_t vcount_before = vertices(composite).size();
      size_t fcount_before = faces(composite).size();
      
      try {
        for(auto v : vertices(meshes[i])) {
          auto new_v = composite.add_vertex(vpm_src[v]);
          vmap[v] = new_v;
        }
        
        for(auto f : faces(meshes[i])) {
          std::vector<typename boost::graph_traits<Mesh>::vertex_descriptor> face_verts;
          for(auto v : vertices_around_face(halfedge(f, meshes[i]), meshes[i])) {
            face_verts.push_back(vmap[v]);
          }
          if(face_verts.size() >= 3) {
            composite.add_face(face_verts);
          }
        }
        
        size_t vcount_after = vertices(composite).size();
        size_t fcount_after = faces(composite).size();
        std::cerr << "  Appended: added " << (vcount_after - vcount_before) << " vertices, " 
                  << (fcount_after - fcount_before) << " faces\n";
      } catch(const std::exception& e) {
        std::cerr << "  Error during append: " << e.what() << "\n";
        // If append fails, just keep temp (skip this mesh)
        composite = temp;
        consecutive_failures++;  // Count append failures too
      }
    }
  }
  
  return composite;
}

void print_usage(const char* program_name) {
  std::cout << "Usage: " << program_name << " <critical_mesh.off> [--cut cut_mesh.off | --fill fill_mesh.off]\n";
  std::cout << "\nAt least one of --cut or --fill must be specified. Both can be provided.\n";
  std::cout << "\nBehavior:\n";
  std::cout << "  - If 1st optional mesh is --cut: Top mode with critical mesh\n";
  std::cout << "  - If 1st optional mesh is --fill: Bottom mode with critical mesh\n";
  std::cout << "  - If 2nd optional mesh is --fill: Bottom mode with each selected critical part\n";
  std::cout << "  - If 2nd optional mesh is --cut: Top mode with each selected critical part\n";
  std::cout << "\nExamples:\n";
  std::cout << "  " << program_name << " critical.off --cut cut.off\n";
  std::cout << "  " << program_name << " critical.off --fill fill.off\n";
  std::cout << "  " << program_name << " critical.off --cut cut.off --fill fill.off\n";
  std::cout << "  " << program_name << " critical.off --fill fill.off --cut cut.off\n";
}

int main(int argc, char* argv[])
{
  std::cout << "=== Lift Layer Mesh Operations ===\n\n";
  
  if(argc < 4) {
    std::cerr << "Error: Insufficient arguments.\n\n";
    print_usage(argv[0]);
    return 1;
  }

  // Parse command line arguments
  std::string critical_file = argv[1];
  std::string first_mesh_file;
  std::string second_mesh_file;
  std::string first_mesh_type;  // "cut" or "fill"
  std::string second_mesh_type; // "cut" or "fill"
  bool has_first_mesh = false;
  bool has_second_mesh = false;

  for(int i = 2; i < argc; ++i) {
    std::string arg = argv[i];
    if(arg == "--cut" || arg == "--fill") {
      if(i + 1 >= argc) {
        std::cerr << "Error: " << arg << " requires a filename argument.\n";
        return 1;
      }
      
      if(!has_first_mesh) {
        first_mesh_type = (arg == "--cut") ? "cut" : "fill";
        first_mesh_file = argv[i + 1];
        has_first_mesh = true;
        ++i; // Skip the filename
      } else if(!has_second_mesh) {
        second_mesh_type = (arg == "--cut") ? "cut" : "fill";
        second_mesh_file = argv[i + 1];
        has_second_mesh = true;
        ++i; // Skip the filename
      } else {
        std::cerr << "Error: Too many mesh arguments provided.\n";
        return 1;
      }
    } else {
      std::cerr << "Error: Unknown argument '" << arg << "'.\n";
      return 1;
    }
  }

  if(!has_first_mesh) {
    std::cerr << "Error: At least one of --cut or --fill must be specified.\n\n";
    print_usage(argv[0]);
    return 1;
  }

  // Load critical mesh
  Mesh critical_mesh;
  if(!PMP::IO::read_polygon_mesh(critical_file, critical_mesh)) {
    std::cerr << "Error: Cannot read critical mesh from " << critical_file << "\n";
    return 1;
  }
  std::cout << "Loaded critical mesh: " << num_vertices(critical_mesh) << " vertices, " 
            << num_faces(critical_mesh) << " faces\n";

  // Load first mesh
  Mesh first_mesh;
  if(!PMP::IO::read_polygon_mesh(first_mesh_file, first_mesh)) {
    std::cerr << "Error: Cannot read " << first_mesh_type << " mesh from " << first_mesh_file << "\n";
    return 1;
  }
  std::cout << "Loaded " << first_mesh_type << " mesh: " << num_vertices(first_mesh) << " vertices, " 
            << num_faces(first_mesh) << " faces\n";

  // === FIRST OPERATION ===
  std::vector<Mesh> selected_parts_critical;
  std::vector<Mesh> selected_parts_first;
  bool first_is_top = (first_mesh_type == "cut");
  std::string first_prefix = first_is_top ? "1st_top" : "1st_bottom";
  
  std::cout << "\n=== First Operation: " << (first_is_top ? "TOP" : "BOTTOM") 
            << " mode with critical mesh and " << first_mesh_type << " mesh ===\n";
  
  PMP::split_and_select_surface_meshes<Mesh, K>(
    critical_mesh, first_mesh, first_is_top,
    selected_parts_critical, selected_parts_first, first_prefix
  );
  
  std::cout << "First operation results:\n";
  std::cout << "  Selected critical parts: " << selected_parts_critical.size() << "\n";
  std::cout << "  Selected " << first_mesh_type << " parts: " << selected_parts_first.size() << "\n";

  // === SECOND OPERATION (if specified) ===
  if(has_second_mesh) {
    // Load second mesh
    Mesh second_mesh;
    if(!PMP::IO::read_polygon_mesh(second_mesh_file, second_mesh)) {
      std::cerr << "Error: Cannot read " << second_mesh_type << " mesh from " << second_mesh_file << "\n";
      return 1;
    }
    std::cout << "\nLoaded " << second_mesh_type << " mesh: " << num_vertices(second_mesh) << " vertices, " 
              << num_faces(second_mesh) << " faces\n";

    bool second_is_top = (second_mesh_type == "cut");
    std::string second_mode = second_is_top ? "TOP" : "BOTTOM";
    
    std::cout << "\n=== Second Operation: " << second_mode 
              << " mode with each selected critical part and " << second_mesh_type << " mesh ===\n";

    // Process each selected critical part from the first operation
    std::vector<Mesh> all_second_selected_critical;
    std::vector<Mesh> all_second_selected_second;
    
    for(size_t i = 0; i < selected_parts_critical.size(); ++i) {
      std::cout << "\nProcessing selected critical part " << i << " of " << selected_parts_critical.size() << "...\n";
      
      std::vector<Mesh> second_selected_critical;
      std::vector<Mesh> second_selected_second;
      std::string second_prefix = "2nd_" + std::string(second_is_top ? "top" : "bottom") + "_" + std::to_string(i);
      
      PMP::split_and_select_surface_meshes<Mesh, K>(
        selected_parts_critical[i], second_mesh, second_is_top,
        second_selected_critical, second_selected_second, second_prefix
      );
      
      std::cout << "  Part " << i << " results: " 
                << second_selected_critical.size() << " critical, "
                << second_selected_second.size() << " " << second_mesh_type << "\n";
      
      // Accumulate results
      all_second_selected_critical.insert(
        all_second_selected_critical.end(),
        second_selected_critical.begin(),
        second_selected_critical.end()
      );
      all_second_selected_second.insert(
        all_second_selected_second.end(),
        second_selected_second.begin(),
        second_selected_second.end()
      );
    }
    
    std::cout << "\nSecond operation total results:\n";
    std::cout << "  Total selected critical parts: " << all_second_selected_critical.size() << "\n";
    std::cout << "  Total selected " << second_mesh_type << " parts: " << all_second_selected_second.size() << "\n";
    
    // Create composite mesh: second_selected_second + second_selected_critical + selected_parts_first
    std::cout << "\n=== Creating Composite Mesh (Second Operation) ===\n";
    std::vector<Mesh> meshes_to_join;
    
    // Add all second operation results
    meshes_to_join.insert(meshes_to_join.end(), 
                          all_second_selected_second.begin(), 
                          all_second_selected_second.end());
    meshes_to_join.insert(meshes_to_join.end(), 
                          all_second_selected_critical.begin(), 
                          all_second_selected_critical.end());
    // Add first operation's first mesh results (exclude selected_parts_critical)
    meshes_to_join.insert(meshes_to_join.end(), 
                          selected_parts_first.begin(), 
                          selected_parts_first.end());
    
    std::cout << "Joining " << meshes_to_join.size() << " mesh parts...\n";
    Mesh composite = join_meshes(meshes_to_join);
    
    std::cout << "Composite mesh: " << num_vertices(composite) << " vertices, " 
              << num_faces(composite) << " faces\n";
    
    CGAL::IO::write_polygon_mesh("composite.off", composite, CGAL::parameters::stream_precision(17));
    std::cout << "Composite mesh written to composite.off\n";
    
  } else {
    // Create composite mesh: selected_parts_critical + selected_parts_first
    std::cout << "\n=== Creating Composite Mesh (First Operation Only) ===\n";
    std::vector<Mesh> meshes_to_join;
    
    meshes_to_join.insert(meshes_to_join.end(), 
                          selected_parts_critical.begin(), 
                          selected_parts_critical.end());
    meshes_to_join.insert(meshes_to_join.end(), 
                          selected_parts_first.begin(), 
                          selected_parts_first.end());
    
    std::cout << "Joining " << meshes_to_join.size() << " mesh parts...\n";
    Mesh composite = join_meshes(meshes_to_join);
    
    std::cout << "Composite mesh: " << num_vertices(composite) << " vertices, " 
              << num_faces(composite) << " faces\n";
    
    CGAL::IO::write_polygon_mesh("composite.off", composite, CGAL::parameters::stream_precision(17));
    std::cout << "Composite mesh written to composite.off\n";
  }

  std::cout << "\n=== Processing Complete ===\n";
  std::cout << "All output files have been written to the current directory.\n";
  
  return 0;
}
