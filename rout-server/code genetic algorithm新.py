"""
@author: Bojun Li
"""

from time import time
# set start time
start_time = time()

from rasterio import open as rio_open
from pathlib import Path
from pyproj import Geod, Transformer, CRS
from pickle import dump, load, HIGHEST_PROTOCOL
from osmnx import graph_from_file
from networkx import set_node_attributes, astar_path
from rtree import index
from geopandas import GeoSeries,read_file
from matplotlib.patches import Patch
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import subplots, savefig, Line2D
from shapely.geometry import LineString
from numpy.random import randint, uniform
from math import sqrt
from heapq import heappush, heappop
from itertools import count
import networkx as nx
from networkx.algorithms.shortest_paths.weighted import _weight_function

def astar_path_1(G, source, target, ww, heuristic=None, weight="weight"):

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
                h = heuristic(neighbor, target,ww)
            enqueued[neighbor] = ncost, h
            push(queue, (ncost + h, next(c), neighbor, ncost, curnode))
            
        traveled_path.append(graph.nodes(data=True)[curnode])
        
        azF, azB, distance_horizontal_path = g.inv(traveled_path[-2]['x'], traveled_path[-2]['y'], traveled_path[-1]['x'], traveled_path[-1]['y'])
        horizontal_path_length.append(distance_horizontal_path)

        vertical_path_distance.append(abs(traveled_path[-2]["Elevation"]-traveled_path[-1]["Elevation"]))


    raise nx.NetworkXNoPath(f"Node {target} not reachable from {source}")

def distance (x1,y1,x2,y2):
    return sqrt((x1-x2)**2+(y1-y2)**2)
    

def crossover(parents, offspring_size):
    """
     * Single point crossover function for decimal numbers
    """

    # loop enough times to make the offspring size
    offspring = []
    for i in range(offspring_size):
        # get binary representations of mum and dad's y values
        
        # get binary representations of mum and dad's y values
        parent_1 = list('{:029b}'.format(parents[(i) % len(parents)]['y']))
        parent_2 = list('{:029b}'.format(parents[(i+1) % len(parents)]['y']))

        # swap some random chromasomes (bits) in the binary strings
        for r in randint(len(parent_1)-1, size=len(parent_1) // 2):
            parent_1[r] = parent_2[r]

        # convert back to number and store in a dictionary
        offspring.append({'y': int("".join(parent_1), 2), 'fitness': None})

    # return the next generation
    return offspring


def mutation(population, mutation_probability, max_mutation):
    

    """
    * Mutate a value by +/- max_mutation
    """
    # mutation changes a single gene in each offspring randomly.
    for i in range(len(population)):

        # does this child want to mutate?
        if (uniform() < mutation_probability):

            # apply the random value as a mutation to the child
            uniform_value = int(uniform(-max_mutation, max_mutation))
            
            
            # if not greater than 0 then store the value of the variation in population[i]['y']
            population[i]['y'] += uniform_value

            # determine if the value after mutation is less than 0
            if  population[i]['y'] < 0:
                population[i]['y'] = population[i]['y'] - uniform_value 
                population[i]['y'] = population[i]['y'] + int(uniform(population[i]['y'], max_mutation))


            elif population[i]['y'] > digital:
                
                 # if it is less than zero, mutate the value of population[i]['y'] from 0 to max_mutation
                population[i]['y'] = population[i]['y'] - uniform_value 

                population[i]['y'] = population[i]['y'] + int(uniform(-max_mutation,digital - population[i]['y']))

    # return the resulting offspring
    return population

def heuristic_function(a, b,ww):

    # extract nodes
    start = graph.nodes(data=True)[a]
    end = graph.nodes(data=True)[b]

    # compute forward and back azimuths, plus distance to gain the distance between cuurent point to the destination
    azF, azB, distance_aD = g.inv(start['x'], start['y'], end['x'], end['y'])

    # compute forward and back azimuths, plus distance to gain the distance between orgin to the cuurent point
    azF, azB, distance_Oa = g.inv(graph.nodes(data=True)[fromNode]['x'], graph.nodes(data=True)[fromNode]['y'], start['x'], start['y'])

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
    vertical_inf = sum(vertical_path_distance)/(abs(graph.nodes(data=True)[fromNode]["Elevation"]-start["Elevation"])+0.01**10)
   
    # calculate the heuristic_score
    heuristic_score =  ((horizontal_inf * ww) + (vertical_inf * (digital-ww)))*distance_aD**1
    

    del traveled_path[-1]
    del horizontal_path_length[-1]
    del vertical_path_distance[-1]

    return heuristic_score


def get_fitness(population):
    
    # loop each individual in the population
    for individual in population:

        # store the argument of individual['y'] to w for later heuristic function calls
        ww = individual['y']
        distance_h = distance_v = 0

        # calculate the shortest avoid sickness path across the network and extract from graph
        get_path = astar_path_1(graph, source = fromNode, target = toNode,ww =ww, heuristic = heuristic_function)
   
            # loop through each node in the shortest path and load into list
        line = []
        for path_node in get_path:
    
            # get the relevant node from the graph with lat lng data
            node = graph.nodes(data=True)[path_node]
    
            # load the lat lng data into the lineString
            line.append([node['x'], node['y']])
        
        for i in range(len(line)-1):
            distance_h = distance(line[i][0],line[i][1],line[i+1][0],line[i+1][1]) + distance_h
            distance_v += abs(graph.nodes(data=True)[list(idx.nearest((line[i][0], line[i][1]), 1))[0]]["Elevation"]-graph.nodes(data=True)[list(idx.nearest((line[i+1][0], line[i+1][1]), 1))[0]]["Elevation"])
        
        h_d_OD = distance(line[0][0],line[0][1],line[len(line)-1][0],line[len(line)-1][1])
        v_d_OD = abs(graph.nodes(data=True)[list(idx.nearest((line[0][0], line[0][1]), 1))[0]]["Elevation"]-graph.nodes(data=True)[list(idx.nearest((line[len(line)-1][0], line[len(line)-1][1]), 1))[0]]["Elevation"])
        s_h = distance_h/h_d_OD
        s_v = distance_v/v_d_OD


        # update fitness as the followd function
        #individual['fitness'] =round( 1/digital*((sum(horizontal_path)/(distance_OD) * individual['y']) + (sum(vertical_path_distance)/(abs(graph.nodes(data=True)[fromNode]["Elevation"]-graph.nodes(data=True)[toNode]["Elevation"])) * (digital-individual['y']))),4)
        # individual['fitness'] = round(1/digital*(s_h* w+s_v*(digital-w)),5)
        individual['fitness'] = 1/digital*(s_h* w+s_v*(digital-w))
        # individual['fitness'] = 1/digital*(s_h+s_v)


        # clear the data of traveled_path list except for the first point(origin) to avoid repeat operations in heuristic_function
        del traveled_path[1:]

        # clear the data of horizontal_path and vertical_path_distance list to avoid repeat operations in heuristic_function
        horizontal_path_length.clear()
        vertical_path_distance.clear()
        get_path.clear()

    return population



''' setting the parameter'''

# settings genetic algorithm parameter
digital                = 100000000       	# number of digits of w
pop_size                = 50              	# population size
num_parents_mating      = 10            	# mating pool size (how many of the pop get to breed)
threshold               = 1.1            	# the desired precision of the result
mutation_probability    = 0.1             	# probability of a child mutating
max_mutation            = 0.1*digital   	# 0.1 digital max mutation
w                       = 0               	# defining weights

        
''' open road graph, building and road shapefiles for map'''

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

# if no building pickle data are available
if not Path('./data/building_shape.pkl').is_file():
    print("No building pickle file found, loading data (takes a long time).")

    # load the data into a Python object from the building file
    building = read_file('./data/building/building.shp')

    # pickle the results for future use
    with open('./data/building_shape.pkl', 'wb') as output:
        dump(building, output, HIGHEST_PROTOCOL)

else:
    print("Building pickle file found successfully, loading data.")

    # extract data from pickle file
    with open('./data/building_shape.pkl', 'rb') as input:
        building = load(input)


# if no road pickle data are available
if not Path('./data/road_shape.pkl').is_file():
    print("No road pickle file found, loading data (takes a long time).")

    # load the data into a Python object from the road file
    road = read_file('./data/road/road.shp')

    # pickle the results for future use
    with open('./data/road_shape.pkl', 'wb') as output:
        dump(road, output, HIGHEST_PROTOCOL)

else:

    print("Road pickle file found successfully, loading data.")

    # extract data from pickle file
    with open('./data/road_shape.pkl', 'rb') as input:
        road = load(input)
        
''' The main code part of the calculation path '''

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

# calculate the 'from' and 'to' node as the nearest to the specified coordinates
# fromNode = list(idx.nearest((-3.1349487320403573, 54.600570793995885), 1))[0]
# toNode = list(idx.nearest((-3.120014192245436, 54.597885857614706), 1))[0]
fromNode = list(idx.nearest((-3.135738372802735,54.597645915375686), 1))[0]
toNode = list(idx.nearest((-3.1217908859252934,54.603838407460245), 1))[0]


# fromNode = list(idx.nearest((-3.128271102905274,54.60018271275991), 1))[0]
# toNode = list(idx.nearest((-3.135008811950684,54.600729844417415), 1))[0]
azF, azB, distance_OD = g.inv(graph.nodes(data=True)[fromNode]['x'], graph.nodes(data=True)[fromNode]['y'], graph.nodes(data=True)[toNode]['x'], graph.nodes(data=True)[toNode]['y'])

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
    
    # creat the empty list to sotre the traveled point
    traveled_path= []

    # creat the empty list to sotre the traveled point dem data
    horizontal_path_length = []

    # creat the empty list to sotre the change in dem value between two points
    vertical_path_distance = []

    # add the origin to the traveled_path list
    traveled_path.append(graph.nodes(data=True)[fromNode])
    
    print('Estimation of parameters using genetic algorithms')
   
    # create the initial population (array of dictionaries) then calculate the fitness for each individual
    population = get_fitness([{'y': int(y), 'fitness': None} for y in uniform(low=0, high=digital, size=pop_size)])
 
    # initialise loop varibles
    generation = 0
    best_fit = float("inf")
    previous_best_fit = None

    # loop until we either find a solution to within the threshold, or the solutions stop improving
    while best_fit > threshold and best_fit != previous_best_fit:

        # select the best parents in the population for mating
        parents = sorted(population, key=lambda individual: individual['fitness'], reverse=False)[:num_parents_mating]

        # get the next generation, mutate and update fitness values
        population = get_fitness(mutation(crossover(parents, pop_size), mutation_probability, max_mutation))
       
        # get the current best individual
        best_match = sorted(population, key=lambda individual: individual['fitness'], reverse=False)[0]
        previous_best_fit = best_fit
        best_fit = best_match['fitness']

        # increment generation counter and report current fitness (1 = perfect)
        generation += 1
        print(f"\tgeneration {generation}: {best_fit}")

    # report the best match and output
    print(f"Best weighting: {best_match['y']} (fitness: {best_match['fitness']} generations: {generation})")
    ww = best_match['y']
  
    
    print('Planning routes to avoid motion sickness')
    
    # calculate the best avoid car sickness path across the network and extract from graph
    get_path = astar_path_1(graph, source = fromNode, target = toNode,ww =ww, heuristic = heuristic_function)

      
    # loop through each node in the shortest path and load into list
    line = []
    for path_node in get_path:

        # get the relevant node from the graph with lat lng data
        node = graph.nodes(data=True)[path_node]

        # load the lat lng data into the lineString
        line.append([node['x'], node['y']])

    # store as a LineString
    lineString = LineString(line)

# define the project of UTM zone 34
utm34 = "+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert linestring to GeoSeries and project to UTM zone 34
path = GeoSeries(lineString, crs="+proj=longlat +datum=WGS84 +no_defs").to_crs(utm34)


# create map axis object, remove axes, set title
fig, my_ax = subplots(1, 1, figsize=(16, 10))
my_ax.axis('off')
my_ax.set(title="The path of avoid car sickness")

# set bounds
buffer = 300
my_ax.set_xlim([path.geometry.iloc[0].bounds[0] - buffer, path.geometry.iloc[0].bounds[2] + buffer])
my_ax.set_ylim([path.geometry.iloc[0].bounds[1] - buffer, path.geometry.iloc[0].bounds[3] + buffer])

#add the road
road.to_crs(utm34).plot(
    ax=my_ax,
    color='#a6cee3',
    linewidth = 2,
    )

# add the buildings
building.to_crs(utm34).plot(
    ax=my_ax,
    color='grey',
    linewidth = 1,
    )

# add the path
path.plot(
    ax=my_ax,
    color='red',
    linewidth = 4,
    )

# manually draw a legend
my_ax.legend([
    Patch(facecolor='grey', label='Buildings'),
    Line2D([0], [0], color='#a6cee3', lw=2),
    Line2D([1], [1], color='red', lw=2)],
    ['Buildings', 'Road', 'Path'], loc='lower left')

# add north arrow
x, y, arrow_length = 0.99, 0.99, 0.1
my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
               arrowprops=dict(facecolor='black', width=5, headwidth=15),
               ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

# add scalebar
my_ax.add_artist(ScaleBar(dx=1, units="m", location="lower right"))

# save the result
savefig(f'out/path.png', bbox_inches='tight')

# report the end
print("done!")

# report runtime
print(f"completed in: {time() - start_time} seconds")	# NO CODE BELOW HERE      
