# coding:utf-8
from flask import Flask,render_template,request
from rasterio import open as rio_open
from pathlib import Path
from pyproj import Geod, Transformer, CRS
from pickle import dump, load, HIGHEST_PROTOCOL
from osmnx import graph_from_file
from networkx import set_node_attributes, astar_path, NetworkXNoPath
from rtree import index
from heapq import heappush, heappop
from itertools import count
import networkx as nx
from networkx.algorithms.shortest_paths.weighted import _weight_function
from math import sqrt


def distance (x1,y1,x2,y2):
    return sqrt((x1-x2)**2+(y1-y2)**2)

def astar_path_1(G, source, target, heuristic=None, weight="weight"):

    if source not in G or target not in G:
        msg = f"Either source {source} or target {target} is not in G"
        raise nx.NodeNotFound(msg)

    if heuristic is None:
        # The default heuristic is h=0 - same as Dijkstra's algorithm
        def heuristic(u, v):
            return 0

    push = heappush
    pop = heappop
    weight = _weight_function(G, weight)

    # The queue stores priority, node, cost to reach, and parent.
    # Uses Python heapq to keep in priority order.
    # Add a counter to the queue to prevent the underlying heap from
    # attempting to compare the nodes themselves. The hash breaks ties in the
    # priority and is guaranteed unique for all nodes in the graph.
    c = count()
    queue = [(0, next(c), source, 0, None)]

    # Maps enqueued nodes to distance of discovered paths and the
    # computed heuristics to target. We avoid computing the heuristics
    # more than once and inserting the node into the queue too many times.
    enqueued = {}
    # Maps explored nodes to parent closest to the source.
    explored = {}

    while queue:
        # Pop the smallest item from queue.
        _, __, curnode, dist, parent = pop(queue)

        if curnode == target:
            path = [curnode]
            node = parent
            while node is not None:
                path.append(node)
                node = explored[node]
            path.reverse()
            return path

        if curnode in explored:
            # Do not override the parent of starting node
            if explored[curnode] is None:
                continue

            # Skip bad paths that were enqueued before finding a better one
            qcost, h = enqueued[curnode]
            if qcost < dist:
                continue

        explored[curnode] = parent
        
        for neighbor, w in G[curnode].items():
            ncost = dist + weight(curnode, neighbor, w)
            if neighbor in enqueued:
                qcost, h = enqueued[neighbor]
                # if qcost <= ncost, a less costly path from the
                # neighbor to the source was already determined.
                # Therefore, we won't attempt to push this neighbor
                # to the queue
                if qcost <= ncost:
                    continue
            else:
                h = heuristic(neighbor, target,source)
            enqueued[neighbor] = ncost, h
            push(queue, (ncost + h, next(c), neighbor, ncost, curnode))
            
        # add the current point to the traveled_path list
        traveled_path.append(graph.nodes(data=True)[curnode])
        
        # calculate the distance between the current point and the previous point,
        azF, azB, distance_horizontal_path = g.inv(traveled_path[-2]['x'], traveled_path[-2]['y'], traveled_path[-1]['x'], traveled_path[-1]['y'])
        
        # add the calculated distance to the horizontal_path_length list
        horizontal_path_length.append(distance_horizontal_path)

        # calculate the absolute value of the elevation difference between the current point and the previous point, and add this value to the vertical_path_distance list
        vertical_path_distance.append(abs(traveled_path[-2]["Elevation"]-traveled_path[-1]["Elevation"]))

    raise nx.NetworkXNoPath(f"Node {target} not reachable from {source}")

def heuristic_function_1(a, b):

    # extract nodes
    start = graph.nodes(data=True)[a]
    end = graph.nodes(data=True)[b]

    # compute forward and back azimuths, plus distance to gain the distance between cuurent point to the destination
    azF, azB, distance_aD = g.inv(start['x'], start['y'], end['x'], end['y'])

    return distance_aD

def heuristic_function_2(a, b,source):

    # extract nodes
    start = graph.nodes(data=True)[a]
    end = graph.nodes(data=True)[b]

    # compute forward and back azimuths, plus distance to gain the distance between cuurent point to the destination
    azF, azB, distance_aD = g.inv(start['x'], start['y'], end['x'], end['y'])

    # compute forward and back azimuths, plus distance to gain the distance between orgin to the cuurent point
    azF, azB, distance_Oa = g.inv(graph.nodes(data=True)[source]['x'], graph.nodes(data=True)[source]['y'], start['x'], start['y'])

    # store the current star point in the traveled_path list
    traveled_path.append(start)
    
    # compute forward and back azimuths, plus distance to gain the distance between the last point to the cuurent point
    azF, azB, distance_horizontal_path = g.inv(traveled_path[-2]['x'], traveled_path[-2]['y'], start['x'], start['y'])

    # store the distance between the current point and the previous point in the horizontal_path list
    horizontal_path_length.append(distance_horizontal_path)

    # store the absolute value of the elevation difference between the current point and the previous point
    vertical_path_distance.append(abs(traveled_path[-2]["Elevation"]-start["Elevation"]))

    # calculate the horizontal influence:
    # divide the distance of the path that has been traversed by the straight-line distance from the starting point to the current point to get the value of sinousity
    horizontal_inf = sum(horizontal_path_length)/(distance_Oa+0.01**10)

    # calculate the vertical influence:
    vertical_inf = sum(vertical_path_distance)/(abs(graph.nodes(data=True)[source]["Elevation"]-start["Elevation"])+0.01**10)
   
    # calculate the heuristic_score
    heuristic_score =  ((horizontal_inf * wh) + (vertical_inf * wv))*distance_aD**30
    del traveled_path[-1]
    del horizontal_path_length[-1]
    del vertical_path_distance[-1]
    
    return heuristic_score

def get_shortest_path(origin,destination):
    
    # calculate the 'from' and 'to' node as the nearest to the specified coordinates
    fromNode = list(idx.nearest(origin, 1))[0]
    toNode = list(idx.nearest(destination, 1))[0]
    
    # use try statement to catch exceptions
    try: 
    
        # calculate the shortest path across the network and extract from graph
        get_path = astar_path(graph, source = fromNode, target = toNode, heuristic = heuristic_function_1)
    # catch exception for no path available and exit the script
    except NetworkXNoPath:
            # print("Sorry, there is no path between those locations in the provided network")
        exit()
    # loop through each node in the shortest path and load into list
    line = []   
    for path_node in get_path:

        # get the relevant node from the graph with lat lng data
        node = graph.nodes(data=True)[path_node]

        # load the lat lng data into the lineString
        line.append([node['y'], node['x']])

    return line

def get_carsickness_avoid_path(origin,destination):
    
    # calculate the 'from' and 'to' node as the nearest to the specified coordinates
    fromNode = list(idx.nearest((origin), 1))[0]
    toNode = list(idx.nearest((destination), 1))[0]

    #Open the dtm data of keswick
    with rio_open ('./data/dem.tif') as d:
        
        # read the data out of band 1 in the dataset
        dem_data = d.read(1)
    
        # create empty dictionaries for store dem data
        dem_dict = {}
    
        # loop through the nodes in the OSM graph
        for graph_node in graph.nodes(data=True):
    
            # transform the coordinates from geo to projected
            x, y = transformer.transform(graph_node[1]['x'], graph_node[1]['y'], direction='FORWARD')
    
            #fill the dictionary to the dem in image space
            dem_dict[graph_node[0]] = dem_data[d.index(x, y)]
    
        # set the node attributes to elevation from the dictionary and call it "Elevation"
        set_node_attributes(graph, dem_dict, "Elevation")

    
        # add the origin to the traveled_path list
        traveled_path.append(graph.nodes(data=True)[fromNode])
        
        # use try statement to catch exceptions
        try: 
    
            # calculate the shortest path across the network and extract from graph
            get_path = astar_path_1(graph, source = fromNode, target = toNode, heuristic = heuristic_function_2)
    
        # catch exception for no path available and exit the script
        except NetworkXNoPath:
            # print("Sorry, there is no path between those locations in the provided network")
            exit()
            
    # loop through each node in the shortest path and load into list
        line = []
        for path_node in get_path:
    
            # get the relevant node from the graph with lat lng data
            node = graph.nodes(data=True)[path_node]
    
            # load the lat lng data into the lineString
            line.append([node['y'], node['x']])
    
    return line

''' open road graph for map'''

# if no osm pickle data are available
if not Path('./data/osm.pkl').is_file():
    print("No osm pickle file found, loading data (takes a long time).")

    # load the data into a Python object from the osm file
    graph = graph_from_file('./data/keswick.osm', name="keswick")

    # pickle the results for future use
    with open('./data/osm.pkl', 'wb') as output:
        dump(graph, output, HIGHEST_PROTOCOL)

else:
    print("OSM pickle file found successfully, loading data.")

    # extract data from pickle file
    with open('./data/osm.pkl', 'rb') as input:
        graph = load(input)

#weighting vhorizontal influence by user define
wh = 0.24221997

#weighting vertical influence by user define
wv = 1-wh

# set which ellipsoid you would like to use
g = Geod(ellps='WGS84')
proj_string  = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs'
geo_string = "+proj=latlong +datum=WGS84 "

# initialise a PyProj Transformer to transform coordinates
transformer = Transformer.from_crs(CRS.from_proj4(geo_string), CRS.from_proj4(proj_string), always_xy=True)

# create spatial index from graph
idx = index.Index()
for id, data in list(graph.nodes(data=True)):
    idx.insert(int(id), (data['x'], data['y'], data['x'], data['y']))
 # creat the empty list to sotre the traveled point
traveled_path= []

# creat the empty list to sotre the traveled point dem data
horizontal_path_length = []

# creat the empty list to sotre the change in dem value between two points
vertical_path_distance = []

# creat the empty list to sotre the start point spatial index value
fromNode = []

# creat the empty list to sotre the end point spatial index value
toNode = []

app = Flask(__name__)

# Binding routes in the form of decorators
@app.route('/')
def hello_world():
     return render_template('index1.html')
 
# Path planning
@app.route('/showRoueData')
def showRoueData():
    
    # get start point coordinates from the web click
    startpoint =eval(request.args.get("start"))
    
    # get the end point coordinates from the web click
    endpoint = eval(request.args.get("end"))

    # define navigation type,1 means speed priority, 2 means avoiding motion sickness
    nav_type = float(request.args.get("type"))
    print("0-------",startpoint,endpoint,nav_type)

    route = []

    if nav_type == 1:

        route = get_shortest_path(startpoint, endpoint)
        
    elif nav_type == 2:

        route = get_carsickness_avoid_path(startpoint, endpoint)
    
    # call the planning method according to the category, and replace the planning result with the data below to achieve the rendering of the route
    return {
        "success":True,
        "message":"Route planning success",
        "data":route    }

if __name__ == '__main__':

    # view routing information in the entire flask through url_map
    print (app.url_map)

    # start the flask program
    app.run()
