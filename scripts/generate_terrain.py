#!/usr/bin/env python3
"""
Generate test terrain meshes with various embankment types.

Usage:
    python generate_terrain.py <type> <count> [output_file] [--grid-size N] [--cell-size S]

Types:
    circular_mound    - Circular mounds/hills
    ridge_horizontal  - Horizontal ridges (along X axis)
    ridge_vertical    - Vertical ridges (along Y axis)
    ridge_diagonal    - Diagonal ridges (45 degrees)
    mixed             - Mix of different embankment types

Options:
    --parallel        Generate evenly-spaced parallel ridges with uniform size

Examples:
    python generate_terrain.py circular_mound 5
    python generate_terrain.py ridge_vertical 3
    python generate_terrain.py ridge_horizontal 4 my_terrain.off
    python generate_terrain.py mixed 6 --grid-size 150
    python generate_terrain.py ridge_horizontal 4 --parallel --grid-size 30
"""

import argparse
import math
import random
import sys
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Embankment:
    """Base class for embankment parameters"""
    cx: float  # center x
    cy: float  # center y
    height: float
    slope_width: float


@dataclass
class CircularMound(Embankment):
    radius: float


@dataclass
class LinearRidge(Embankment):
    angle: float  # radians
    half_length: float
    top_width: float


def embankment_profile(dist_from_center: float, width: float, height: float, slope_width: float) -> float:
    """Create a trapezoidal embankment profile."""
    if dist_from_center < width / 2.0:
        return height
    elif dist_from_center < width / 2.0 + slope_width:
        t = (dist_from_center - width / 2.0) / slope_width
        return height * (1.0 - t)
    return 0.0


def circular_mound_height(x: float, y: float, mound: CircularMound) -> float:
    """Calculate height contribution from a circular mound."""
    dist = math.sqrt((x - mound.cx) ** 2 + (y - mound.cy) ** 2)
    return embankment_profile(dist, mound.radius * 2.0, mound.height, mound.slope_width)


def linear_ridge_height(x: float, y: float, ridge: LinearRidge) -> float:
    """Calculate height contribution from a linear ridge."""
    dx = x - ridge.cx
    dy = y - ridge.cy
    
    # Rotate to local coordinates
    local_along = dx * math.cos(ridge.angle) + dy * math.sin(ridge.angle)
    local_across = -dx * math.sin(ridge.angle) + dy * math.cos(ridge.angle)
    
    # Check if within length
    if abs(local_along) > ridge.half_length + ridge.slope_width:
        return 0.0
    
    # Taper at ends
    length_factor = 1.0
    if abs(local_along) > ridge.half_length:
        length_factor = 1.0 - (abs(local_along) - ridge.half_length) / ridge.slope_width
    
    return embankment_profile(abs(local_across), ridge.top_width, ridge.height, ridge.slope_width) * length_factor


def generate_random_mounds(count: int, terrain_size: float, margin: float = 15.0) -> List[CircularMound]:
    """Generate random circular mounds."""
    mounds = []
    for _ in range(count):
        mounds.append(CircularMound(
            cx=random.uniform(margin, terrain_size - margin),
            cy=random.uniform(margin, terrain_size - margin),
            height=random.uniform(4.0, 12.0),
            slope_width=random.uniform(4.0, 8.0),
            radius=random.uniform(4.0, 10.0)
        ))
    return mounds


def generate_random_ridges(count: int, terrain_size: float, angle: float, margin: float = 20.0) -> List[LinearRidge]:
    """Generate random linear ridges with specified angle."""
    ridges = []
    for _ in range(count):
        ridges.append(LinearRidge(
            cx=random.uniform(margin, terrain_size - margin),
            cy=random.uniform(margin, terrain_size - margin),
            height=random.uniform(5.0, 10.0),
            slope_width=random.uniform(4.0, 7.0),
            angle=angle,
            half_length=random.uniform(15.0, 35.0),
            top_width=random.uniform(3.0, 6.0)
        ))
    return ridges


def generate_parallel_ridges(
    count: int, 
    terrain_size: float, 
    angle: float,
    height: float = 6.0,
    top_width: float = 2.0,
    slope_width: float = 3.0,
    spacing: float = None  # explicit spacing between ridges (center to center)
) -> List[LinearRidge]:
    """Generate evenly-spaced parallel ridges with uniform size.
    
    If spacing is provided, ridges are placed with that exact distance between centers.
    Otherwise, ridges are evenly distributed across the terrain.
    """
    ridges = []
    
    # Ridge half-length spans almost the full terrain
    margin = terrain_size * 0.1
    half_length = (terrain_size - 2 * margin) / 2.0 * 0.9
    
    # Calculate positions
    if count == 1:
        spacing_positions = [terrain_size / 2.0]
    elif spacing is not None:
        # Use explicit spacing, center the group
        total_span = spacing * (count - 1)
        start_pos = (terrain_size - total_span) / 2.0
        spacing_positions = [start_pos + spacing * i for i in range(count)]
    else:
        # Evenly distribute across terrain
        margin_size = terrain_size * 0.1
        usable_size = terrain_size - 2 * margin_size
        step = usable_size / (count + 1)
        spacing_positions = [margin_size + step * (i + 1) for i in range(count)]
    
    for i, pos in enumerate(spacing_positions):
        # Determine center based on ridge orientation
        if abs(angle) < 0.01:  # Horizontal ridge
            cx = terrain_size / 2.0
            cy = pos
        elif abs(angle - math.pi / 2.0) < 0.01:  # Vertical ridge
            cx = pos
            cy = terrain_size / 2.0
        else:  # Diagonal
            cx = terrain_size / 2.0
            cy = pos
        
        ridges.append(LinearRidge(
            cx=cx,
            cy=cy,
            height=height,
            slope_width=slope_width,
            angle=angle,
            half_length=half_length,
            top_width=top_width
        ))
    
    return ridges


def generate_mixed(count: int, terrain_size: float, margin: float = 20.0) -> Tuple[List[CircularMound], List[LinearRidge]]:
    """Generate a mix of mounds and ridges."""
    mounds = []
    ridges = []
    
    for i in range(count):
        embankment_type = random.choice(['mound', 'ridge_h', 'ridge_v', 'ridge_d'])
        
        if embankment_type == 'mound':
            mounds.append(CircularMound(
                cx=random.uniform(margin, terrain_size - margin),
                cy=random.uniform(margin, terrain_size - margin),
                height=random.uniform(4.0, 12.0),
                slope_width=random.uniform(4.0, 8.0),
                radius=random.uniform(4.0, 10.0)
            ))
        else:
            if embankment_type == 'ridge_h':
                angle = 0.0
            elif embankment_type == 'ridge_v':
                angle = math.pi / 2.0
            else:
                angle = random.choice([math.pi / 4.0, -math.pi / 4.0])
            
            ridges.append(LinearRidge(
                cx=random.uniform(margin, terrain_size - margin),
                cy=random.uniform(margin, terrain_size - margin),
                height=random.uniform(5.0, 10.0),
                slope_width=random.uniform(4.0, 7.0),
                angle=angle,
                half_length=random.uniform(15.0, 35.0),
                top_width=random.uniform(3.0, 6.0)
            ))
    
    return mounds, ridges


def generate_terrain_mesh(
    grid_size: int,
    cell_size: float,
    mounds: List[CircularMound],
    ridges: List[LinearRidge]
) -> Tuple[List[Tuple[float, float, float]], List[Tuple[int, int, int]]]:
    """Generate terrain mesh vertices and faces."""
    
    # Generate height field
    heights = []
    for i in range(grid_size + 1):
        row = []
        for j in range(grid_size + 1):
            x = i * cell_size
            y = j * cell_size
            z = 0.0
            
            # Add contributions from all mounds
            for mound in mounds:
                z = max(z, circular_mound_height(x, y, mound))
            
            # Add contributions from all ridges
            for ridge in ridges:
                z = max(z, linear_ridge_height(x, y, ridge))
            
            row.append(z)
        heights.append(row)
    
    # Create vertices
    vertices = []
    for i in range(grid_size + 1):
        for j in range(grid_size + 1):
            x = i * cell_size
            y = j * cell_size
            z = heights[i][j]
            vertices.append((x, y, z))
    
    # Create triangular faces
    def vertex_index(i, j):
        return i * (grid_size + 1) + j
    
    faces = []
    for i in range(grid_size):
        for j in range(grid_size):
            v00 = vertex_index(i, j)
            v10 = vertex_index(i + 1, j)
            v01 = vertex_index(i, j + 1)
            v11 = vertex_index(i + 1, j + 1)
            
            faces.append((v00, v10, v11))
            faces.append((v00, v11, v01))
    
    return vertices, faces


def write_off_file(filename: str, vertices: List[Tuple[float, float, float]], faces: List[Tuple[int, int, int]]):
    """Write mesh to OFF file format."""
    with open(filename, 'w') as f:
        f.write("OFF\n")
        f.write(f"{len(vertices)} {len(faces)} 0\n")
        
        for v in vertices:
            f.write(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        
        for face in faces:
            f.write(f"3 {face[0]} {face[1]} {face[2]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate test terrain meshes with various embankment types.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python generate_terrain.py circular_mound 5
    python generate_terrain.py ridge_vertical 3
    python generate_terrain.py ridge_horizontal 4 my_terrain.off
    python generate_terrain.py mixed 6 --grid-size 150
    python generate_terrain.py ridge_horizontal 4 --parallel --grid-size 30
        """
    )
    
    parser.add_argument('type', choices=[
        'circular_mound', 'ridge_horizontal', 'ridge_vertical', 'ridge_diagonal', 'mixed'
    ], help='Type of embankment to generate')
    
    parser.add_argument('count', type=int, help='Number of embankments to generate')
    
    parser.add_argument('output', nargs='?', default=None,
                        help='Output OFF file (default: terrain_<type>_<count>.off)')
    
    parser.add_argument('--grid-size', type=int, default=100,
                        help='Grid size (default: 100)')
    
    parser.add_argument('--cell-size', type=float, default=1.0,
                        help='Cell size (default: 1.0)')
    
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility')
    
    parser.add_argument('--parallel', action='store_true',
                        help='Generate evenly-spaced parallel ridges with uniform size')
    
    parser.add_argument('--height', type=float, default=6.0,
                        help='Ridge/mound height for parallel mode (default: 6.0)')
    
    parser.add_argument('--top-width', type=float, default=2.0,
                        help='Ridge top width for parallel mode (default: 2.0)')
    
    parser.add_argument('--slope-width', type=float, default=3.0,
                        help='Slope width for parallel mode (default: 3.0)')
    
    parser.add_argument('--spacing', type=float, default=None,
                        help='Spacing between parallel ridges (center to center). If not set, ridges are evenly distributed.')
    
    args = parser.parse_args()
    
    if args.seed is not None:
        random.seed(args.seed)
    
    terrain_size = args.grid_size * args.cell_size
    
    mounds = []
    ridges = []
    
    if args.type == 'circular_mound':
        mounds = generate_random_mounds(args.count, terrain_size)
        type_desc = f"{args.count} circular mound(s)"
    
    elif args.type == 'ridge_horizontal':
        if args.parallel:
            ridges = generate_parallel_ridges(
                args.count, terrain_size, angle=0.0,
                height=args.height, top_width=args.top_width, slope_width=args.slope_width,
                spacing=args.spacing
            )
            type_desc = f"{args.count} parallel horizontal ridge(s)"
        else:
            ridges = generate_random_ridges(args.count, terrain_size, angle=0.0)
            type_desc = f"{args.count} horizontal ridge(s)"
    
    elif args.type == 'ridge_vertical':
        if args.parallel:
            ridges = generate_parallel_ridges(
                args.count, terrain_size, angle=math.pi / 2.0,
                height=args.height, top_width=args.top_width, slope_width=args.slope_width,
                spacing=args.spacing
            )
            type_desc = f"{args.count} parallel vertical ridge(s)"
        else:
            ridges = generate_random_ridges(args.count, terrain_size, angle=math.pi / 2.0)
            type_desc = f"{args.count} vertical ridge(s)"
    
    elif args.type == 'ridge_diagonal':
        if args.parallel:
            ridges = generate_parallel_ridges(
                args.count, terrain_size, angle=math.pi / 4.0,
                height=args.height, top_width=args.top_width, slope_width=args.slope_width,
                spacing=args.spacing
            )
            type_desc = f"{args.count} parallel diagonal ridge(s)"
        else:
            ridges = generate_random_ridges(args.count, terrain_size, angle=math.pi / 4.0)
            type_desc = f"{args.count} diagonal ridge(s)"
    
    elif args.type == 'mixed':
        mounds, ridges = generate_mixed(args.count, terrain_size)
        type_desc = f"mixed terrain with {len(mounds)} mound(s) and {len(ridges)} ridge(s)"
    
    # Generate mesh
    vertices, faces = generate_terrain_mesh(args.grid_size, args.cell_size, mounds, ridges)
    
    # Determine output filename
    if args.output:
        output_file = args.output
    else:
        output_file = f"terrain_{args.type}_{args.count}.off"
    
    # Write to file
    write_off_file(output_file, vertices, faces)
    
    # Print summary
    print(f"Generated terrain: {type_desc}")
    print(f"  Grid size: {args.grid_size} x {args.grid_size}")
    print(f"  Terrain size: {terrain_size} x {terrain_size}")
    print(f"  Vertices: {len(vertices)}")
    print(f"  Faces: {len(faces)}")
    print()
    
    if mounds:
        print("Circular mounds:")
        for i, m in enumerate(mounds, 1):
            print(f"  {i}. center=({m.cx:.1f}, {m.cy:.1f}), radius={m.radius:.1f}, height={m.height:.1f}")
    
    if ridges:
        print("Linear ridges:")
        for i, r in enumerate(ridges, 1):
            angle_deg = math.degrees(r.angle)
            print(f"  {i}. center=({r.cx:.1f}, {r.cy:.1f}), angle={angle_deg:.0f}°, length={r.half_length*2:.1f}, height={r.height:.1f}")
    
    print(f"\nSaved to: {output_file}")


if __name__ == "__main__":
    main()
