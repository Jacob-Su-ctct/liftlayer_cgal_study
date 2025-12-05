#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef K::Vector_3 Vector_3;
typedef K::Ray_3 Ray_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

// Compute base elevation (minimum Z or average of lowest points)
double compute_base_elevation(const Surface_mesh& mesh, double percentile = 0.1)
{
    std::vector<double> z_values;
    for (auto v : mesh.vertices()) {
        z_values.push_back(mesh.point(v).z());
    }
    std::sort(z_values.begin(), z_values.end());
    
    size_t idx = static_cast<size_t>(z_values.size() * percentile);
    return z_values[idx];
}

// Compute mesh centroid (horizontal center)
std::pair<double, double> compute_horizontal_center(const Surface_mesh& mesh)
{
    double sum_x = 0, sum_y = 0;
    size_t count = 0;
    
    for (auto v : mesh.vertices()) {
        const Point_3& p = mesh.point(v);
        sum_x += p.x();
        sum_y += p.y();
        count++;
    }
    
    return {sum_x / count, sum_y / count};
}

// Structure to hold vertex expansion info
struct VertexExpansionInfo {
    Point_3 original_pos;
    Vector_3 normal;
    double horiz_nx, horiz_ny;  // normalized horizontal direction
    double horiz_len;           // horizontal component magnitude
    bool is_slope;
    double max_allowed_expansion;
};

// Find the nearest opposing slope point to limit expansion
double find_max_expansion(const Point_3& pos, 
                          double dir_x, double dir_y,
                          const std::vector<std::pair<Point_2, std::pair<double, double>>>& slope_points,
                          double default_max)
{
    double min_dist = default_max;
    
    for (const auto& sp : slope_points) {
        const Point_2& other_pos = sp.first;
        double other_dir_x = sp.second.first;
        double other_dir_y = sp.second.second;
        
        // Vector from this point to other point (horizontal)
        double dx = other_pos.x() - pos.x();
        double dy = other_pos.y() - pos.y();
        double dist = std::sqrt(dx * dx + dy * dy);
        
        if (dist < 1e-6) continue;
        
        // Normalize
        double ndx = dx / dist;
        double ndy = dy / dist;
        
        // Check if other point is in our expansion direction
        double dot_our_dir = dir_x * ndx + dir_y * ndy;
        
        // Check if other point's expansion direction is toward us (opposing)
        double dot_their_dir = other_dir_x * (-ndx) + other_dir_y * (-ndy);
        
        // If both slopes face each other (opposing directions)
        if (dot_our_dir > 0.3 && dot_their_dir > 0.3) {
            // Limit expansion to half the distance (meet in the middle)
            double allowed = (dist / 2.0) * 0.9;  // 0.9 for safety margin
            if (allowed < min_dist) {
                min_dist = allowed;
            }
        }
    }
    
    return min_dist;
}

// Expand embankments with collision detection
void expand_embankments_with_collision_detection(Surface_mesh& mesh,
                                                   double expansion_factor,
                                                   double base_z)
{
    namespace PMP = CGAL::Polygon_mesh_processing;
    
    // Store original positions and compute info
    std::map<Surface_mesh::Vertex_index, VertexExpansionInfo> vertex_info;
    
    for (auto v : mesh.vertices()) {
        VertexExpansionInfo info;
        info.original_pos = mesh.point(v);
        info.is_slope = false;
        info.max_allowed_expansion = expansion_factor;
        vertex_info[v] = info;
    }
    
    // Compute vertex normals
    auto vnormals = mesh.add_property_map<Surface_mesh::Vertex_index, Vector_3>(
        "v:normal", Vector_3(0, 0, 1)).first;
    PMP::compute_vertex_normals(mesh, vnormals);
    
    const double slope_threshold = 0.1;
    
    // First pass: identify slopes and their directions
    std::vector<std::pair<Point_2, std::pair<double, double>>> slope_points;
    
    for (auto v : mesh.vertices()) {
        Vector_3 n = vnormals[v];
        double nx = n.x();
        double ny = n.y();
        double horiz_len = std::sqrt(nx * nx + ny * ny);
        
        vertex_info[v].normal = n;
        vertex_info[v].horiz_len = horiz_len;
        
        if (horiz_len >= slope_threshold) {
            vertex_info[v].is_slope = true;
            vertex_info[v].horiz_nx = nx / horiz_len;
            vertex_info[v].horiz_ny = ny / horiz_len;
            
            const Point_3& p = vertex_info[v].original_pos;
            slope_points.push_back({
                Point_2(p.x(), p.y()),
                {vertex_info[v].horiz_nx, vertex_info[v].horiz_ny}
            });
        }
    }
    
    std::cout << "Found " << slope_points.size() << " slope points\n";
    
    // Second pass: compute max allowed expansion for each slope vertex
    for (auto v : mesh.vertices()) {
        if (!vertex_info[v].is_slope) continue;
        
        double max_exp = find_max_expansion(
            vertex_info[v].original_pos,
            vertex_info[v].horiz_nx,
            vertex_info[v].horiz_ny,
            slope_points,
            expansion_factor
        );
        
        vertex_info[v].max_allowed_expansion = max_exp;
    }
    
    // Third pass: apply expansion with limits
    int limited_count = 0;
    for (auto v : mesh.vertices()) {
        if (!vertex_info[v].is_slope) continue;
        
        const Point_3& p = vertex_info[v].original_pos;
        double expansion = std::min(expansion_factor, vertex_info[v].max_allowed_expansion);
        
        if (expansion < expansion_factor) {
            limited_count++;
        }
        
        double new_x = p.x() + vertex_info[v].horiz_nx * expansion;
        double new_y = p.y() + vertex_info[v].horiz_ny * expansion;
        
        mesh.point(v) = Point_3(new_x, new_y, p.z());
    }
    
    std::cout << "Limited expansion for " << limited_count << " vertices to avoid overlap\n";
    
    mesh.remove_property_map(vnormals);
}

// Compute intersection of two 2D lines (slope expansion lines)
// Line 1: point p1 with direction d1
// Line 2: point p2 with direction d2
// Returns true if intersection found, sets int_x, int_y
bool compute_line_intersection_2d(
    double p1_x, double p1_y, double d1_x, double d1_y,
    double p2_x, double p2_y, double d2_x, double d2_y,
    double& int_x, double& int_y)
{
    // Solve: p1 + t*d1 = p2 + s*d2
    // t*d1_x - s*d2_x = p2_x - p1_x
    // t*d1_y - s*d2_y = p2_y - p1_y
    
    double denom = d1_x * (-d2_y) - d1_y * (-d2_x);
    if (std::abs(denom) < 1e-10) {
        return false;  // Parallel lines
    }
    
    double dx = p2_x - p1_x;
    double dy = p2_y - p1_y;
    
    double t = (dx * (-d2_y) - dy * (-d2_x)) / denom;
    
    int_x = p1_x + t * d1_x;
    int_y = p1_y + t * d1_y;
    
    return true;
}

// Alternative: Fill valleys instead of limiting expansion
// This expands slopes fully and handles crossing slopes by computing proper intersections
void expand_embankments_fill_valleys(Surface_mesh& mesh,
                                      double expansion_factor,
                                      double base_z)
{
    namespace PMP = CGAL::Polygon_mesh_processing;
    
    // Store original positions
    std::map<Surface_mesh::Vertex_index, Point_3> original_positions;
    for (auto v : mesh.vertices()) {
        original_positions[v] = mesh.point(v);
    }
    
    // Compute vertex normals
    auto vnormals = mesh.add_property_map<Surface_mesh::Vertex_index, Vector_3>(
        "v:normal", Vector_3(0, 0, 1)).first;
    PMP::compute_vertex_normals(mesh, vnormals);
    
    const double slope_threshold = 0.1;
    
    // Collect slope info for each vertex
    struct SlopeInfo {
        double horiz_nx, horiz_ny;
        double horiz_len;
        bool is_slope;
        bool is_flat_top;
        bool is_flat_base;
        Point_3 new_pos;
        bool position_updated;
    };
    std::map<Surface_mesh::Vertex_index, SlopeInfo> slope_info;
    
    // Find max height for reference
    double max_height = 0;
    for (auto v : mesh.vertices()) {
        double h = original_positions[v].z() - base_z;
        if (h > max_height) max_height = h;
    }
    
    // First pass: identify slopes and compute expansion directions
    std::cout << "Pass 1: Identifying slopes...\n";
    for (auto v : mesh.vertices()) {
        Vector_3 n = vnormals[v];
        double horiz_len = std::sqrt(n.x() * n.x() + n.y() * n.y());
        double z = original_positions[v].z();
        
        SlopeInfo si;
        si.horiz_len = horiz_len;
        si.is_slope = (horiz_len >= slope_threshold);
        si.is_flat_top = (!si.is_slope && z > base_z + max_height * 0.8);
        si.is_flat_base = (!si.is_slope && z <= base_z + 0.1);
        si.new_pos = original_positions[v];
        si.position_updated = false;
        
        if (si.is_slope) {
            si.horiz_nx = n.x() / horiz_len;
            si.horiz_ny = n.y() / horiz_len;
        } else {
            si.horiz_nx = 0;
            si.horiz_ny = 0;
        }
        slope_info[v] = si;
    }
    
    int slope_count = 0, flat_top_count = 0, flat_base_count = 0;
    for (auto& kv : slope_info) {
        if (kv.second.is_slope) slope_count++;
        if (kv.second.is_flat_top) flat_top_count++;
        if (kv.second.is_flat_base) flat_base_count++;
    }
    std::cout << "  Slopes: " << slope_count << ", Flat tops: " << flat_top_count 
              << ", Flat base: " << flat_base_count << "\n";
    
    // Second pass: compute initial expanded positions for slope vertices
    std::cout << "Pass 2: Computing expanded positions...\n";
    for (auto v : mesh.vertices()) {
        if (!slope_info[v].is_slope) continue;
        
        const Point_3& p = original_positions[v];
        double nx = slope_info[v].horiz_nx;
        double ny = slope_info[v].horiz_ny;
        
        // Compute expanded position
        double new_x = p.x() + nx * expansion_factor;
        double new_y = p.y() + ny * expansion_factor;
        slope_info[v].new_pos = Point_3(new_x, new_y, p.z());
    }
    
    // Third pass: detect crossing slopes and compute intersection points
    // For each pair of opposing slopes that would cross, compute where their
    // expansion lines intersect and use that as the meeting point
    std::cout << "Pass 3: Detecting crossings and computing intersections...\n";
    
    std::set<Surface_mesh::Vertex_index> crossing_vertices;
    std::map<Surface_mesh::Vertex_index, Point_3> intersection_positions;
    
    for (auto v : mesh.vertices()) {
        if (!slope_info[v].is_slope) continue;
        
        const Point_3& p = original_positions[v];
        double nx = slope_info[v].horiz_nx;
        double ny = slope_info[v].horiz_ny;
        
        // Check for crossing with opposing slopes
        for (auto v2 : mesh.vertices()) {
            if (v >= v2) continue;  // Avoid duplicate checks
            if (!slope_info[v2].is_slope) continue;
            
            const Point_3& p2 = original_positions[v2];
            double nx2 = slope_info[v2].horiz_nx;
            double ny2 = slope_info[v2].horiz_ny;
            
            // Check if directions are opposing
            double dot = nx * nx2 + ny * ny2;
            if (dot > -0.5) continue;  // Not opposing
            
            // Check if they are at similar Z levels (part of same valley)
            if (std::abs(p.z() - p2.z()) > max_height * 0.5) continue;
            
            // Check if the expanded positions would cross
            const Point_3& new_p = slope_info[v].new_pos;
            const Point_3& new_p2 = slope_info[v2].new_pos;
            
            double dist_before = std::sqrt((p.x()-p2.x())*(p.x()-p2.x()) + 
                                           (p.y()-p2.y())*(p.y()-p2.y()));
            double dist_after = std::sqrt((new_p.x()-new_p2.x())*(new_p.x()-new_p2.x()) + 
                                          (new_p.y()-new_p2.y())*(new_p.y()-new_p2.y()));
            
            // If they would cross or get very close
            if (dist_after < dist_before) {
                crossing_vertices.insert(v);
                crossing_vertices.insert(v2);
                
                // Compute intersection of the two expansion lines
                double int_x, int_y;
                if (compute_line_intersection_2d(
                        p.x(), p.y(), nx, ny,
                        p2.x(), p2.y(), nx2, ny2,
                        int_x, int_y)) {
                    
                    // Check if intersection is in the expansion direction for both
                    double t1 = (int_x - p.x()) * nx + (int_y - p.y()) * ny;
                    double t2 = (int_x - p2.x()) * nx2 + (int_y - p2.y()) * ny2;
                    
                    if (t1 > 0 && t2 > 0) {
                        // Valid intersection - both moving toward it
                        // Limit to expansion_factor distance
                        double dist1 = std::sqrt((int_x - p.x())*(int_x - p.x()) + 
                                                 (int_y - p.y())*(int_y - p.y()));
                        double dist2 = std::sqrt((int_x - p2.x())*(int_x - p2.x()) + 
                                                 (int_y - p2.y())*(int_y - p2.y()));
                        
                        // Move to intersection or max expansion, whichever is less
                        if (dist1 <= expansion_factor) {
                            intersection_positions[v] = Point_3(int_x, int_y, p.z());
                        }
                        if (dist2 <= expansion_factor) {
                            intersection_positions[v2] = Point_3(int_x, int_y, p2.z());
                        }
                    }
                }
            }
        }
    }
    
    std::cout << "  Found " << crossing_vertices.size() << " crossing vertices\n";
    
    // Apply intersection positions to crossing vertices
    for (auto& kv : intersection_positions) {
        slope_info[kv.first].new_pos = kv.second;
        slope_info[kv.first].position_updated = true;
    }
    
    // Fourth pass: identify faces to remove
    // Remove faces where all vertices have moved to approximately the same position
    // (these are the valley floor faces that get "crushed" by opposing slopes)
    std::cout << "Pass 4: Identifying faces to remove...\n";
    
    std::set<Surface_mesh::Face_index> faces_to_remove;
    
    for (auto f : mesh.faces()) {
        std::vector<Surface_mesh::Vertex_index> face_verts;
        for (auto v : mesh.vertices_around_face(mesh.halfedge(f))) {
            face_verts.push_back(v);
        }
        
        if (face_verts.size() != 3) continue;
        
        // Check if this face is in the valley (all vertices at base level originally)
        bool all_at_base = true;
        for (auto v : face_verts) {
            if (original_positions[v].z() > base_z + max_height * 0.3) {
                all_at_base = false;
                break;
            }
        }
        
        if (!all_at_base) continue;
        
        // Check if any vertex is a crossing vertex
        bool has_crossing = false;
        for (auto v : face_verts) {
            if (crossing_vertices.find(v) != crossing_vertices.end()) {
                has_crossing = true;
                break;
            }
        }
        
        if (!has_crossing) continue;
        
        // Check if the face would become degenerate (vertices too close together)
        const Point_3& np0 = slope_info[face_verts[0]].new_pos;
        const Point_3& np1 = slope_info[face_verts[1]].new_pos;
        const Point_3& np2 = slope_info[face_verts[2]].new_pos;
        
        // Compute new face area (2D projection)
        double ax = np1.x() - np0.x(), ay = np1.y() - np0.y();
        double bx = np2.x() - np0.x(), by = np2.y() - np0.y();
        double new_area = std::abs(ax * by - ay * bx) / 2.0;
        
        // Original area
        const Point_3& op0 = original_positions[face_verts[0]];
        const Point_3& op1 = original_positions[face_verts[1]];
        const Point_3& op2 = original_positions[face_verts[2]];
        double oax = op1.x() - op0.x(), oay = op1.y() - op0.y();
        double obx = op2.x() - op0.x(), oby = op2.y() - op0.y();
        double orig_area = std::abs(oax * oby - oay * obx) / 2.0;
        
        // If face shrinks significantly or inverts, remove it
        if (orig_area > 1e-6 && new_area < orig_area * 0.1) {
            faces_to_remove.insert(f);
        }
    }
    
    std::cout << "  Marked " << faces_to_remove.size() << " faces for removal\n";
    
    // Fifth pass: apply new positions to all slope vertices
    std::cout << "Pass 5: Applying new positions...\n";
    for (auto v : mesh.vertices()) {
        if (slope_info[v].is_slope) {
            mesh.point(v) = slope_info[v].new_pos;
        }
    }
    
    // Sixth pass: remove valley faces
    std::cout << "Pass 6: Removing valley faces...\n";
    for (auto f : faces_to_remove) {
        if (mesh.is_valid(f) && !mesh.is_removed(f)) {
            CGAL::Euler::remove_face(mesh.halfedge(f), mesh);
        }
    }
    
    // Collect garbage to clean up removed elements
    mesh.collect_garbage();
    
    // Seventh pass: fill holes created by face removal
    std::cout << "Pass 7: Filling holes...\n";
    std::vector<Surface_mesh::Halfedge_index> border_cycles;
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    
    std::cout << "  Found " << border_cycles.size() << " boundary cycles (holes)\n";
    
    for (auto he : border_cycles) {
        // Count cycle size
        int cycle_size = 0;
        auto curr = he;
        do {
            cycle_size++;
            curr = mesh.next(curr);
        } while (curr != he && cycle_size < 1000);
        
        // Only fill reasonably-sized holes (not the outer boundary)
        if (cycle_size >= 3 && cycle_size < 100) {
            try {
                PMP::triangulate_hole(mesh, he);
                std::cout << "  Filled hole of size " << cycle_size << "\n";
            } catch (...) {
                std::cout << "  Could not fill hole of size " << cycle_size << "\n";
            }
        }
    }
    
    mesh.remove_property_map(vnormals);
    
    std::cout << "Valley fill complete.\n";
}

int main(int argc, char* argv[])
{
    Surface_mesh mesh;
    
    const char* input_file = (argc > 1) ? argv[1] : "input.off";
    double expansion_factor = (argc > 2) ? std::atof(argv[2]) : 5.0;
    std::string method = (argc > 3) ? argv[3] : "collision";
    
    if (!CGAL::IO::read_polygon_mesh(input_file, mesh)) {
        std::cerr << "Error: Cannot read " << input_file << std::endl;
        return 1;
    }
    
    std::cout << "Loaded mesh: " << mesh.number_of_vertices() 
              << " vertices, " << mesh.number_of_faces() << " faces\n";
    
    double base_z = compute_base_elevation(mesh, 0.1);
    std::cout << "Base elevation: " << base_z << "\n";
    
    double max_z = base_z;
    for (auto v : mesh.vertices()) {
        if (mesh.point(v).z() > max_z) {
            max_z = mesh.point(v).z();
        }
    }
    double max_height = max_z - base_z;
    std::cout << "Max height above base: " << max_height << "\n";
    
    if (method == "collision") {
        std::cout << "Using collision detection method\n";
        expand_embankments_with_collision_detection(mesh, expansion_factor, base_z);
    } else if (method == "valley") {
        std::cout << "Using valley fill method\n";
        expand_embankments_fill_valleys(mesh, expansion_factor, base_z);
    } else {
        std::cerr << "Unknown method: " << method << "\n";
        std::cerr << "Use 'collision' or 'valley'\n";
        return 1;
    }
    
    std::cout << "Embankments expanded by factor " << expansion_factor << "\n";
    
    CGAL::IO::write_polygon_mesh("output_expanded.off", mesh);
    std::cout << "Saved to output_expanded.off\n";
    
    return 0;
}