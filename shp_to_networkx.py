import geopandas as gpd
import shapely as sh
import networkx as nx
from collections import defaultdict
import pyproj
from functools import partial
import matplotlib.pyplot as plt

file = "./test-geometry/test-geometry.shp"

# projections for coordinates (lat/lon) and distances
# https://epsg.io/4326
crs_gps = "+proj=longlat +datum=WGS84 +no_defs"

# https://epsg.io/3654
crs_distance = "+proj=tmerc +lat_0=41.08333333333334 +lon_0=-71.5 +k=0.99999375 +x_0=99999.99998983997 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=us-ft +no_defs"

buffersize = 5e-11  # set the buffer around the line to search for intersections.  should be < 1e-6

# get GeoDataFrame of linestrings
gdf = gpd.read_file(file)
gdf = gdf.to_crs(crs_gps)
# initialize the main graph
G = nx.Graph()

# helper structures for converting intersections and endpoints to nodes
node_counter = 0
pt_to_node = {}


def node_tuple(pt):
    """return the cached node_id and rounded off point for a given coord"""
    precision = 8
    global pt_to_node
    global node_counter
    pt_tuple = (round(pt[0], precision),
                round(pt[1], precision))
    if pt_tuple not in pt_to_node:
        pt_to_node[pt_tuple] = node_counter
        node_counter += 1
    return pt_to_node[pt_tuple], pt_tuple

# helper structure for tying node ids to rows in the geo dataframe
row_to_nodes = defaultdict(set)


def add_to_graph(coord, row, row2=None):
    """centralize logic to add a new node to the graph 
    and tie it to a row in the dataframe"""
    global G
    global node_to_row
    node_id, node_pt = node_tuple(coord)
    G.add_node(node_id, point=sh.geometry.Point(node_pt))
    row_to_nodes[row].add(node_id)
    if row2 is not None:
        row_to_nodes[row2].add(node_id)

for row1 in range(len(gdf)):
    cur_line = gdf.iloc[row1].geometry
    cur_line_buffered = cur_line.buffer(buffersize)
    # add the start and end of this linestring as nodes
    cur_start = cur_line.coords[0]
    cur_end = cur_line.coords[-1]
    add_to_graph(cur_start, row1)
    add_to_graph(cur_end, row1)
    for row2 in range(row1+1, len(gdf)):
        other_line = gdf.iloc[row2]
        # check for intersection or near-touching
        if other_line.geometry.intersects(cur_line_buffered):
            overlap = other_line.geometry.intersection(cur_line_buffered)
            try:
                overlap_list = list(overlap)
            except TypeError:  # only had a single intersection
                overlap_list = [overlap]
            overlap_pts = [obj.centroid.coords[0] for obj in overlap_list]
            for pt in overlap_pts:
                add_to_graph(pt, row1, row2)


# project the coordinates to get useful distances
# http://toblerity.org/shapely/manual.html#other-transformations

transform_point = partial(pyproj.transform,
                          pyproj.Proj(crs_gps),
                          pyproj.Proj(crs_distance))

# now the graph contains nodes for all endpoints and intersections...
# add edges across the nodes which comprise a linesegment from the original gdf
for row, nodes in row_to_nodes.items():
    cur_line = sh.ops.transform(transform_point,
                                gdf.iloc[row].geometry)
    projected_distances = [
        (id, cur_line.project(sh.ops.transform(transform_point,
                                               G.nodes[id]['point']),
                              normalized=False))
        for id in nodes]
    projected_distances.sort(key=lambda x: x[1])
    for ix in range(1, len(projected_distances)):
        node1, dist1 = projected_distances[ix-1]
        node2, dist2 = projected_distances[ix-2]
        G.add_edge(node1, node2,
                   distance=round(dist2-dist1),
                   row=row,
                   linestring=gdf.iloc[row].geometry)

# plot for debugging
for u, v in G.edges:
    u_x, u_y = G.nodes[u]['point'].coords[0]
    v_x, v_y = G.nodes[v]['point'].coords[0]
    plt.plot([u_x, v_x], [u_y, v_y], alpha=.8, linestyle=(0, (1, 10)))

plt.show()
