#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <iostream>
#include <vector>
#include <map>
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

// Alternative: Fill valleys instead of limiting expansion
// This creates a flat bottom where slopes would meet
void expand_embankments_fill_valleys(Surface_mesh& mesh,
                                      double expansion_factor,
                                      double base_z)
{
    namespace PMP = CGAL::Polygon_mesh_processing;
    
    std::map<Surface_mesh::Vertex_index, Point_3> original_positions;
    for (auto v : mesh.vertices()) {
        original_positions[v] = mesh.point(v);
    }
    
    auto vnormals = mesh.add_property_map<Surface_mesh::Vertex_index, Vector_3>(
        "v:normal", Vector_3(0, 0, 1)).first;
    PMP::compute_vertex_normals(mesh, vnormals);
    
    const double slope_threshold = 0.1;
    
    // Collect slope info
    struct SlopeInfo {
        double horiz_nx, horiz_ny;
        bool is_slope;
    };
    std::map<Surface_mesh::Vertex_index, SlopeInfo> slope_info;
    
    for (auto v : mesh.vertices()) {
        Vector_3 n = vnormals[v];
        double horiz_len = std::sqrt(n.x() * n.x() + n.y() * n.y());
        
        SlopeInfo si;
        si.is_slope = (horiz_len >= slope_threshold);
        if (si.is_slope) {
            si.horiz_nx = n.x() / horiz_len;
            si.horiz_ny = n.y() / horiz_len;
        }
        slope_info[v] = si;
    }
    
    // For each vertex, check if expansion would cause overlap
    for (auto v : mesh.vertices()) {
        if (!slope_info[v].is_slope) continue;
        
        const Point_3& p = original_positions[v];
        double nx = slope_info[v].horiz_nx;
        double ny = slope_info[v].horiz_ny;
        
        // Proposed new position
        double new_x = p.x() + nx * expansion_factor;
        double new_y = p.y() + ny * expansion_factor;
        
        // Check for collision with other expanded points
        bool collision = false;
        Surface_mesh::Vertex_index collision_vertex;
        
        for (auto v2 : mesh.vertices()) {
            if (v == v2 || !slope_info[v2].is_slope) continue;
            
            const Point_3& p2 = original_positions[v2];
            
            // Check if directions are opposing
            double dot = nx * slope_info[v2].horiz_nx + ny * slope_info[v2].horiz_ny;
            if (dot > -0.5) continue;  // Not opposing
            
            // Proposed position of v2
            double new_x2 = p2.x() + slope_info[v2].horiz_nx * expansion_factor;
            double new_y2 = p2.y() + slope_info[v2].horiz_ny * expansion_factor;
            
            // Check if they would cross
            double dist_before = std::sqrt((p.x()-p2.x())*(p.x()-p2.x()) + 
                                           (p.y()-p2.y())*(p.y()-p2.y()));
            double dist_after = std::sqrt((new_x-new_x2)*(new_x-new_x2) + 
                                          (new_y-new_y2)*(new_y-new_y2));
            
            if (dist_after < dist_before * 0.3) {  // Getting too close
                collision = true;
                collision_vertex = v2;
                break;
            }
        }
        
        if (collision) {
            // Meet in the middle
            const Point_3& p2 = original_positions[collision_vertex];
            double mid_x = (p.x() + p2.x()) / 2.0;
            double mid_y = (p.y() + p2.y()) / 2.0;
            mesh.point(v) = Point_3(mid_x, mid_y, p.z());
        } else {
            mesh.point(v) = Point_3(new_x, new_y, p.z());
        }
    }
    
    mesh.remove_property_map(vnormals);
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