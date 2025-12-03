#!/usr/bin/env python3
"""
Multiple OFF file viewer using trimesh
Usage: python3 view_multiple_off.py <filename1.off> <filename2.off> ...
"""

import random
import sys
import trimesh
import matplotlib.pyplot as plt
import numpy as np

def view_multiple_off_files(filenames):
    """Load and display multiple OFF files in the same 3D plot"""
    if not filenames:
        print("No files provided.")
        return

    # Create a 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Color palette for different meshes
    colors = ['green', 'lightcoral', 'lightblue', 'lightgray', 'lightyellow', 'lightcyan', 'lightpink', 'lightgreen', 'lavender']
    random_seed = random.randint(0, 9)
    all_vertices = []

    for i, filename in enumerate(filenames):
        try:
            # Load the mesh
            mesh = trimesh.load(filename, process=False)
            all_vertices.append(mesh.vertices)

            # Plot the mesh with a different color
            color = colors[(i+random_seed) % len(colors)]
            if len(mesh.faces) > 0:
                ax.plot_trisurf(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.vertices[:, 2],
                               triangles=mesh.faces, color=color, alpha=0.8, label=filename)
            else:
                # Plot as lines connecting points if no faces (assuming point cloud represents a curve)
                # Sort points by X coordinate to connect them in order
                sorted_indices = np.argsort(mesh.vertices[:, 0])
                sorted_vertices = mesh.vertices[sorted_indices]
                ax.plot(sorted_vertices[:, 0], sorted_vertices[:, 1], sorted_vertices[:, 2],
                       color='red', label=filename, linewidth=2)

        except Exception as e:
            print(f"Error loading {filename}: {e}")
            continue

    if not all_vertices:
        print("No valid meshes loaded.")
        return

    # Combine all vertices to calculate global ranges
    all_vertices = np.vstack(all_vertices)
    x_range = all_vertices[:, 0].max() - all_vertices[:, 0].min()
    y_range = all_vertices[:, 1].max() - all_vertices[:, 1].min()
    z_range = all_vertices[:, 2].max() - all_vertices[:, 2].min()

    # Find the maximum range
    max_range = max(x_range, y_range, z_range)

    # Calculate center points
    x_center = (all_vertices[:, 0].max() + all_vertices[:, 0].min()) / 2
    y_center = (all_vertices[:, 1].max() + all_vertices[:, 1].min()) / 2
    z_center = (all_vertices[:, 2].max() + all_vertices[:, 2].min()) / 2

    # Set equal axis limits to ensure 1:1:1 aspect ratio
    half_range = max_range / 2
    ax.set_xlim(x_center - half_range, x_center + half_range)
    ax.set_ylim(y_center - half_range, y_center + half_range)
    ax.set_zlim(z_center - half_range, z_center + half_range)

    # Ensure the plot is displayed with equal aspect ratio
    ax.set_box_aspect([1, 1, 1])

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Multiple OFF Files')

    # Add legend
    ax.legend()

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 view_multiple_off.py <filename1.off> <filename2.off> ...")
        print("Available OFF files in current directory:")
        import os
        off_files = [f for f in os.listdir('.') if f.endswith('.off')]
        for f in off_files:
            print(f"  {f}")
        sys.exit(1)

    filenames = sys.argv[1:]
    view_multiple_off_files(filenames)