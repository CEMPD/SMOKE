
"""
griddesc2shp.py — Convert SMOKE/CMAQ GRIDDESC to a GeoPackage layer

Overview
    - Parses a SMOKE/CMAQ GRIDDESC file and writes either:
            • A full grid layer (one polygon per cell with ROWNUM/COLNUM)
            • A single extent polygon of the domain
    - Cell vertex coordinates are generated in the native Lambert Conformal Conic (LCC) projection
        and then transformed to WGS84 lon/lat for the output geometry.

Inputs
    - GRIDDESC file containing two sections:
            1) Coordinate systems (terminated by a line containing: ' '  !  end coords.)
            2) Grid definitions that reference a coordinate system
    - Example terminator line in GRIDDESC (must exist):
                ' '  !  end coords.

Outputs
    - A GeoPackage (.gpkg) written to the path you provide (or default). The layer name defaults to
        the grid id, but can be overridden.
    - Attribute fields (full-grid mode):
            • name   – grid id
            • ROWNUM – grid row number (1-based)
            • COLNUM – grid column number (1-based)
            • center_x_m, center_y_m – cell center (meters, native LCC) if --no-centers is NOT used
            • center_lon, center_lat – cell center (WGS84 lon/lat) if --no-centers is NOT used

Projection/Ellipsoid
    - By default this script uses a spherical Earth (a=b=6,370,000 m) for the LCC projection.
        Toggle this via the USE_SPHERICAL_EARTH variable near the top of the file if you need WGS84.

Requirements
    - Python 3.8+
    - Packages: geopandas, shapely, pyproj

Usage (tcsh)
    - Full grid with default output location and layer name:
            python3 griddesc2shp.py -g griddesc_lambertonly_18jan2019_v7.txt -i 12LISTOS

    - Full grid with explicit output path and custom layer name:
            python3 griddesc2shp.py -g /path/to/GRIDDESC.txt -i 12US1 \
                -o /tmp/12US1_domain.gpkg --layer-name 12US1

    - Single extent polygon (no per-cell output):
            python3 griddesc2shp.py -g /path/to/GRIDDESC.txt -i 12US1 --extent-only

    - Skip center attributes in full-grid mode:
            python3 griddesc2shp.py -g /path/to/GRIDDESC.txt -i 12US1 --no-centers

CLI Options
    Required
    - -i/--grid-id           Grid (domain) name/ID to extract (required)

    Optional
    - -g/--griddesc          Path to the GRIDDESC file (default: griddesc_lambertonly_18jan2019_v7.txt)
    - -o/--output            Output GeoPackage path (default: <GRID-ID>_domain.gpkg in current directory)
    - --layer-name           Layer name inside the GeoPackage (default: same as grid id)
    - --extent-only          Write a single extent polygon instead of the full grid (flag; default: off)
    - --no-centers           Skip center attributes in full-grid mode (flag; default: off)

Notes
    - The script expects GRIDDESC numeric fields that may contain Fortran 'D' exponents (e.g., 1.0D+00);
        these are converted automatically.
    - The script prints the native LCC projection string and derived domain corners for quick QA.
"""

import re
import geopandas as gpd
from shapely.geometry import Polygon
import pyproj
import os
import argparse

# Toggle: set to True if you want spherical earth (a=b=6370000 m) instead of WGS84
USE_SPHERICAL_EARTH = True

def _clean_name(raw: str) -> str:
    raw = raw.split('!')[0]  # drop inline comment
    raw = raw.strip()
    raw = re.sub(r"['`,]", "", raw)  # drop quotes/commas/backticks
    raw = re.sub(r"\s+", " ", raw)
    return raw.strip()

def parse_griddesc_all(path: str):
    with open(path, 'r') as f:
        lines = f.readlines()
    try:
        sep_idx = next(i for i, ln in enumerate(lines) if "' '  !  end coords." in ln)
    except StopIteration:
        raise ValueError("Missing coordinate block terminator (' '  !  end coords.)")

    coords = {}
    i = 0
    while i < sep_idx:
        line = lines[i].strip()
        if line.startswith("'") and not line.startswith("!"):
            name = _clean_name(line)
            i += 1
            if i < sep_idx:
                params_line = lines[i].strip()
                nums = [float(tok.replace('D', 'E')) for tok in re.split(r',\s*|\s+', params_line) if tok.strip()]
                coords[name] = nums
        i += 1

    grids = {}
    i = sep_idx + 1
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("'") and not line.startswith("!") and line != "' '":
            gname = _clean_name(line)
            i += 1
            if i < len(lines):
                params_line = lines[i].strip()
                parts = [p.strip() for p in params_line.split(',') if p.strip()]
                if parts:
                    coord_ref = _clean_name(parts[0])
                    rest = [float(p.replace('D', 'E')) for p in parts[1:]]
                    grids[gname] = [coord_ref] + rest
        i += 1
    return coords, grids

def extract_grid(path: str, grid_id: str):
    coords, grids = parse_griddesc_all(path)
    gid_clean = _clean_name(grid_id)
    if gid_clean not in grids:
        # Provide friendly error
        raise ValueError(f"Grid '{grid_id}' not found. Available: {', '.join(sorted(grids.keys()))}")
    grid_params = grids[gid_clean]
    coord_name = grid_params[0]
    if coord_name not in coords:
        raise ValueError(f"Projection '{coord_name}' referenced by grid '{gid_clean}' not defined in coords section.")
    return coords[coord_name], grid_params

def parse_griddesc(file_path, grid_name):
    """
    Parses a GRIDDESC file to extract projection and grid parameters.
    """
    coords = {}
    grids = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Find where coordinate definitions end and grid definitions begin
    try:
        end_coords_index = next(i for i, line in enumerate(lines) if "' '  !  end coords." in line)
    except StopIteration:
        raise ValueError("Could not find the end of the coordinate definitions in the GRIDDESC file.")

    # Parse coordinate definitions
    i = 0
    while i < end_coords_index:
        line = lines[i].strip()
        if line.startswith("'") and not line.startswith("!"):
            coord_name = line.strip("'")
            i += 1
            # Ensure we don't read past the end of the coordinate block
            if i < end_coords_index:
                params_line = lines[i].strip()
                # Clean up Fortran 'D' notation and split
                params = [float(p.replace('D', 'E')) for p in re.split(r',\s*|\s+', params_line)]
                coords[coord_name] = params
        i += 1

    # Parse grid definitions
    i = end_coords_index + 1
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("'") and not line.startswith("!") and line != "' '":
            current_grid_name = line.split('!')[0].strip().strip("'")
            i += 1
            if i < len(lines):
                params_line = lines[i].strip()
                
                # Split by comma, then clean up each part
                parts = [p.strip() for p in params_line.split(',')]
                
                coord_name_part = parts[0].strip("'")
                
                # Handle numeric parts, converting Fortran 'D' notation
                numeric_parts = [float(p.replace('D', 'E')) for p in parts[1:]]
                
                grids[current_grid_name] = [coord_name_part] + numeric_parts
        i += 1
        
    if grid_name not in grids:
        raise ValueError(f"Grid '{grid_name}' not found in GRIDDESC file.")
        
    grid_params = grids[grid_name]
    coord_name = grid_params[0]
    
    if coord_name not in coords:
        raise ValueError(f"Coordinate system '{coord_name}' for grid '{grid_name}' not found.")
        
    coord_params = coords[coord_name]
    
    return coord_params, grid_params

def create_domain_shapefile(
    griddesc_path: str,
    grid_name: str,
    output_gpkg: str,
    full_grid: bool = True,
    add_centers: bool = True,
):
    """Create a GeoPackage layer of a model domain from a GRIDDESC file.

WRITE_NATIVE_LAYER = True    # Also write a layer in native LCC projection for alignment QA
    Simple domain polygon writer.

    full_grid=True  -> one feature per cell (row/col attributes)
    full_grid=False -> single extent polygon.
    add_centers adds center lon/lat & projected meters when full_grid.
    """
    # 1. Parse the file
    #coord_params, grid_params = parse_griddesc(griddesc_path, grid_name)
    coord_params, grid_params = extract_grid(griddesc_path, grid_name)
    
    # Extract projection parameters for LAM_40N97W
    # type, P-alpha, P-beta, P-gamma, xcent, ycent
    proj_type, p_alpha, p_beta, p_gamma, x_cent, y_cent = coord_params
    
    # Extract grid parameters for 12LISTOS
    # coord_name, xorig, yorig, xcell, ycell, ncols, nrows, nthik
    _, xorig, yorig, xcell, ycell, ncols, nrows, _ = grid_params

    # 2. Define the Lambert Conformal projection
    if USE_SPHERICAL_EARTH:
        proj_str = (
            f"+proj=lcc +lat_1={p_alpha} +lat_2={p_beta} +lat_0={y_cent} "
            f"+lon_0={x_cent} +a=6370000.0 +b=6370000.0 +x_0=0 +y_0=0 +units=m +no_defs"
        )
    else:
        proj_str = (
            f"+proj=lcc +lat_1={p_alpha} +lat_2={p_beta} +lat_0={y_cent} "
            f"+lon_0={x_cent} +ellps=WGS84 +datum=WGS84 +x_0=0 +y_0=0 +units=m +no_defs"
        )
    lcc_proj = pyproj.Proj(proj_str)
    # Minimal debug (can be commented out if undesired)
    print("Projection:", proj_str)
    lcc_crs = pyproj.CRS.from_proj4(proj_str)
    print(f"Origin: XORIG={xorig}, YORIG={yorig}; Cell={xcell}x{ycell}; Size={int(ncols)}x{int(nrows)}")
    
    # Define the target projection (WGS84 geographic)
    wgs84_proj = pyproj.Proj(proj='latlong', datum='WGS84')

    # 3. Calculate domain corners in projected coordinates (meters)
    ll_x, ll_y = xorig, yorig
    ur_x, ur_y = xorig + (ncols * xcell), yorig + (nrows * ycell)
    
    corners_proj = {
        'lower_left': (ll_x, ll_y),
        'lower_right': (ur_x, ll_y),
        'upper_right': (ur_x, ur_y),
        'upper_left': (ll_x, ur_y)
    }
    
    print(f"Projected corners (Lambert Conformal, meters) for {grid_name}:")
    for name, (x, y) in corners_proj.items():
        print(f"  {name}: X={x:.2f}, Y={y:.2f}")

    # 4. Transform corners to geographic coordinates (lon/lat)
    transformer = pyproj.Transformer.from_proj(lcc_proj, wgs84_proj, always_xy=True)
    
    corners_latlon = {}
    for name, (x, y) in corners_proj.items():
        lon, lat = transformer.transform(x, y)
        corners_latlon[name] = (lon, lat)

    print("\nTransformed corners (WGS84, lon/lat):")
    for name, (lon, lat) in corners_latlon.items():
        print(f"  {name}: Lon={lon:.6f}, Lat={lat:.6f}")

    if not full_grid:
        # Single extent polygon (lat/lon corners already computed)
        polygon_geom = Polygon([
            corners_latlon['lower_left'],
            corners_latlon['lower_right'],
            corners_latlon['upper_right'],
            corners_latlon['upper_left'],
            corners_latlon['lower_left']
        ])
        gdf = gpd.GeoDataFrame({'name': [grid_name]}, geometry=[polygon_geom], crs='EPSG:4326')
    else:
        features = []
        rows_attr = []
        cols_attr = []
        centers_x = []
        centers_y = []
        centers_lon = []
        centers_lat = []
        for r in range(1, int(nrows)+1):
            y0 = yorig + (r-1) * ycell
            y1 = y0 + ycell
            for c in range(1, int(ncols)+1):
                x0 = xorig + (c-1) * xcell
                x1 = x0 + xcell
                # Build polygon in projected space then transform vertices
                pts_proj = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
                pts_ll = [transformer.transform(px, py) for (px, py) in pts_proj]
                poly = Polygon(pts_ll + [pts_ll[0]])
                features.append(poly)
                rows_attr.append(r)
                cols_attr.append(c)
                if add_centers:
                    cx = x0 + 0.5 * xcell
                    cy = y0 + 0.5 * ycell
                    clon, clat = transformer.transform(cx, cy)
                    centers_x.append(cx)
                    centers_y.append(cy)
                    centers_lon.append(clon)
                    centers_lat.append(clat)
        gdf = gpd.GeoDataFrame(
            {
                'name': [grid_name]*len(features),
                'ROWNUM': rows_attr,
                'COLNUM': cols_attr
            },
            geometry=features,
            crs='EPSG:4326'
        )
        if add_centers:
            gdf['center_x_m'] = centers_x
            gdf['center_y_m'] = centers_y
            gdf['center_lon'] = centers_lon
            gdf['center_lat'] = centers_lat

    gdf.to_file(output_gpkg, driver="GPKG", layer=grid_name)
    msg = f"{len(gdf)} cell polygons" if full_grid else "single extent polygon"
    print(f"Written {msg} to: {output_gpkg}")

def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Create a GeoPackage polygon layer for a SMOKE grid from a GRIDDESC file."
    )
    parser.add_argument(
        '-g', '--griddesc',
        required=True,
        help='Path to GRIDDESC file (e.g.: griddesc.txt)'
    )
    parser.add_argument(
        '-i', '--grid-id',
        required=True,
        help='Grid (domain) name / ID to extract (e.g., 12LISTOS)'
    )
    parser.add_argument(
        '-o', '--output',
        required=False,
        help='Output GeoPackage file path. If omitted: <GRID-ID>_domain.gpkg in current directory.'
    )
    parser.add_argument(
        '--layer-name',
        required=False,
        help='Optional layer name inside GeoPackage (default: same as grid id).'
    )
    parser.add_argument('--extent-only', action='store_true', help='Output only a single extent polygon (default is full grid).')
    parser.add_argument('--no-centers', action='store_true', help='Skip center coordinate attributes (full grid mode).')
    return parser


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    griddesc_path = os.path.abspath(args.griddesc)
    grid_id = args.grid_id
    layer_name = args.layer_name or grid_id

    if args.output:
        output_path = os.path.abspath(args.output)
        # Ensure .gpkg extension unless user explicitly uses another
        if not output_path.lower().endswith('.gpkg'):
            output_path += '.gpkg'
    else:
        output_path = os.path.join(os.getcwd(), f"{grid_id}_domain.gpkg")

    if not os.path.exists(griddesc_path):
        print(f"Error: GRIDDESC file not found at '{griddesc_path}'")
        return 2

    try:
        create_domain_shapefile(
            griddesc_path,
            grid_id,
            output_path,
            full_grid=not args.extent_only,
            add_centers=not args.no_centers,
        )
        # If user wanted a different internal layer name than grid id, rename layer by rewriting.
        if layer_name != grid_id:
            # Read back and re-write with new layer name (GeoPandas cannot rename in place easily)
            gdf = gpd.read_file(output_path, layer=grid_id)
            gdf.to_file(output_path, layer=layer_name, driver='GPKG')
        print(f"Output written: {output_path} (layer: {layer_name})")
    except (ValueError, KeyError) as e:
        print(f"An error occurred: {e}")
        return 1

    return 0


if __name__ == '__main__':
    raise SystemExit(main())

