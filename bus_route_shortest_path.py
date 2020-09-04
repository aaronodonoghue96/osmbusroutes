# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:31:59 2020

@author: aaron

Base map and data from OpenStreetMap and OpenStreetMap Foundation

The shortest path between each pair of nodes and all overall shortest paths produced
are calculated by using the "length" attribute of each edge, which represents the
length of the corresponding road in real life, including any curves on the road.
The shortest path lengths are then calculated from the sum of the length of each 
path between consecutive node pairs.

To reduce the graph size, and ultimately the time and space required to compute the graph,
the coordinates file was examined to determine which areas of Cork contained most of the
bus stops. Once this had been done, the northern, southern, eastern and western extent
of each group of outliers was calculated, and the graphs were created based on these figures,
so that all of the area covered by bus networks would be covered while minimizing wasted
space and minimizing the amount of graphs required.

"""

import tkinter as tk
import osmnx as ox
import osmnx.distance as oxdist
import osmnx.utils as oxutils
import numpy as np
import networkx as nx
import folium as f
import pandas as pd
from math import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import random
import time
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as tiles
import io
from urllib.request import urlopen, Request
from PIL import Image

top = tk.Tk()

WALKING_SPEED = 6 # average walking speed in km/hr
BUS_SPEED = 22 # average bus speed in Cork in km/hr

# a list of bus stops which have charging stations
charging_stations = []

# create a graph of size 12.5 sq km using OSM data from Cork
G1 = ox.graph_from_address("Cork, Cork, Ireland", dist=12500, network_type="drive_service", simplify=False)

# create a graph from all roads within the longitudes of 51.7 and 51.95 degrees north
# and the latitudes of 8.46 and 8.522 degrees west.
# This is to cover the 226 route that leads to Kinsale while minimizing the amount
# of land without bus routes covered
G2 = ox.graph_from_bbox(51.95, 51.7, -8.46, -8.522, network_type="drive_service", simplify=False)

# create a graph of size 5 sq km using OSM data from Ovens
# used to cover Ovens and Ballincollig on the 220 bus route
# while minimizing the amount of land without bus routes covered
G3 = ox.graph_from_address("Ovens, Cork, Ireland", dist=5000, network_type="drive_service", simplify=False)

# combine these three graphs together to cover the areas of Cork with bus routes
# while minimizing the amount of land with no bus routes
G = nx.compose_all([G1, G2, G3])

# a list of roads that must be avoided as they are not on bus routes
# and will otherwise be used by Dijkstra's algorithm as part of shortest paths
roads_to_avoid = []

# remove all possible ways to access the forbidden road from each of its neighbours
for node in G.neighbors(358657):
    roads_to_avoid.append(G[358657][node]) #Robert St (a one-way, narrow street that can't be traversed by bus)

for node in G.neighbors(16353226):
    roads_to_avoid.append(G[16353226][node]) #Jack Lynch Tunnel (neither buses nor pedestrians can travel through here)

# a road in Kinsale that is part of the 226 bus route but has not been added to OSM
# without this road, Dijkstra's algorithm will produce an error
G.add_edge(3963879063, 1901594533, length=26.99)

for road in roads_to_avoid:
    road[0]["length"] = 1000000 #sets the distance to 1 million metres, so it won't be on the shortest path

# a dictionary of all road segments that are part of bus routes
# keys are tuples containing 2 tuples each, indicating the coordinates of the first and last nodes of each road
# values are a list of all bus routes that the road segment indicated by the key is part of
road_segments = {}

# a list of all bus stops used to find the nearest bus stops to the user's specified start and end points
all_bus_stops = []

# used to get OSM map tiles to create the initial selection map for the GUI
# this method comes from markersportal.com
def image_spoof(self, tile): # used to get image of Cork for initial selection
    url = self._image_url(tile) # get the url of the street map API
    req = Request(url) # start request
    req.add_header('User-agent','Anaconda 3') # add user agent to request
    fh = urlopen(req)
    im_data = io.BytesIO(fh.read()) # get image
    fh.close() # close url
    img = Image.open(im_data) # open image with PIL
    img = img.convert(self.desired_tile_form) # set image format
    return img, self.tileextent(tile), 'lower' # reformat for cartopy


print("Graph generated")

# create a map of Cork which will eventually be used to plot the routes,
# bus stops, charging stations and start and end points on
# prefer_canvas=True is used to speed up marker adding time and scrolling in HTML document
map_of_cork = f.Map(prefer_canvas=True)

# used to make each entry in the legend, indicating route name and colour
lgd_txt = '<span style="color: {col};">{txt}</span>'

# dictionary of feature groups
fg_dict = {}

# a timetable of bus routes and bus stops in Cork is used
# to use another city's transport network, replace the string below
# with an Excel file for the city of your choice
timetable = 'timetable-cork.xlsx'

# a dictionary indicating the sheet in the Excel file to be used for each
# bus route's stops, and the colour of each bus route on the map
route_dict = {"226": {"stops": "226-distance", "line_colour": "blue"},
"208": {"stops": "208-distance", "line_colour": "lightblue"},
"207": {"stops": "207-distance", "line_colour": "red"},
"205": {"stops": "205-distance", "line_colour": "green"},
"203": {"stops": "203-distance", "line_colour": "yellow"},
"223": {"stops": "223-distance", "line_colour": "purple"},
"221": {"stops": "221-distance", "line_colour": "orange"},
"220": {"stops": "220-distance", "line_colour": "cyan"},
"216": {"stops": "216-distance", "line_colour": "brown"},
"215": {"stops": "215-distance", "line_colour": "magenta"},
"214": {"stops": "214-distance", "line_colour": "teal"}
}

# give each route a feature group with its number and colour for the legend
# and store in dictionary for later use, i.e. plotting lines and adding bus stops
for key in route_dict.keys():
    fg = f.FeatureGroup(name = lgd_txt.format(txt="Route " + key, col=route_dict[key]["line_colour"]))
    fg_dict[key] = fg

# convert string coordinates to lists, and add each coordinate pair to the list of stops
# which will be used for plotting bus stops in the correct locations
def add_coords(coords, stops):
    for coord in coords:
        pair = [float(s) for s in coord.strip().split(",")]
        stops.append(pair)

# add markers for bus stops on each route to the route's feature group
# so that the bus stops on a route only appear if the route is selected
def add_markers(stops, route_num, charging_stations):

    for i in range(len(stops)):
        # charging stations are found at the first stop (index 0), and every tenth stop after that
        if i % 10 == 0:
            charge = f.Marker(location=[stops[i][0], stops[i][1]], popup=route_num + ' (Charge Point)', icon=f.Icon(color="red"))
            charging_stations += [get_nearest_bus_stop(G, stops[i][0], stops[i][1])]
            fg_dict[route_num].add_child(charge)
        else:
            bus_stop = f.Marker(location=[stops[i][0], stops[i][1]])
            fg_dict[route_num].add_child(bus_stop)

# find the nearest bus stop to the point with longitude y_coord and latitude x_coord
def get_nearest_bus_stop(G, y_coord, x_coord): #adapted from OSMNX's get_nearest_node function
    if len(G) < 1:
        raise ValueError("G must contain at least one node")

    # dump graph node coordinates into a pandas dataframe indexed by node id
    # with x and y columns
    coords = ((p[0], p[1][1], p[1][0]) for p in all_bus_stops) #bus stop number, y coord, x coord
    df = pd.DataFrame(coords, columns=["node", "x", "y"]).set_index("node")

    # add columns to df for the (constant) coordinates of reference point
    df["ref_y"] = y_coord
    df["ref_x"] = x_coord

    # calculate the distance between each node and the reference point
    dists = oxdist.great_circle_vec(lat1=df["ref_y"], lng1=df["ref_x"], lat2=df["y"], lng2=df["x"])

    # nearest node's ID is the index label of the minimum distance
    nearest_node = dists.idxmin()
    return nearest_node
    #find the nearest bus stop to these coords and return that node

#display the shortest path between the user's chosen points using buses where possible
def display_route(route):
    # take a bus where possible, stay on current bus where possible, if going from walk to multiple buses, take first bus
    # initialize bus-related values
    battery = 100
    passengers = 0
    capacity = 50

    # initially, there is no previous or current route, and that will only change once the first bus stop is reached
    prev_route = "N/A"
    curr_route = "N/A"
    
    # create arrays to store X and Y coordinates
    xs = []
    ys = []
    
    # create a path of coordinate pairs
    path = []

    # add each point's coordinates to the arrays
    for p in route:
        xs += [p[1]]
        ys += [p[0]]
        path += [(p[0], p[1])]
    # use a ball tree for fast haversine search on graphs
    nearest = ox.get_nearest_edges(G, xs, ys, 'balltree')

    for i in range(len(route) - 1):
        #only do this is current route is a bus route

        if (route[i], route[i+1]) in road_segments: #if this segment of road is on a bus route
            #if entering a new bus route, transfer battery and passengers over
            prev_route = curr_route # prev_route updated to the current (soon to be previous) route, as curr_route is about to be updated
            curr_route = road_segments[(route[i], route[i+1])][0] # update curr_route to the next route on the path

            if passengers < capacity: # if there is room for more passengers
                on = random.randint(1, capacity - passengers) # add a random amount without exceeding capacity
            else:
                print("Sorry! Bus full!")
                on = 0

            if passengers > 0: # if the bus is not empty
                off = random.randint(1, passengers) # at least one passenger will leave, but cannot have less than 0 passengers
            else:
                print("Bus empty")
                off = 0

            # add the new passengers and remove the leaving passengers from the bus
            passengers += on
            print("%s new passengers have arrived" % on)

            passengers -= off
            print("%s passengers have left" % off)
            print("There are currently %s passengers" % passengers)

            # if a charging station is passed, charge the battery to 100%
            # Future Work: add a time penalty for charging, find the optimal amount to charge,
            # and find a better way to detect charging station bus stops
            if get_nearest_bus_stop(G, route[i][0], route[i][1]) in charging_stations:
                print("Charging battery")
                battery = 100

        # find the current route
        x = G[nearest[i][0]][nearest[i][1]]
        if "name" in x[0].keys(): # print its name if it has one
            print("Travelling on " + str(x[0]["name"]))
        else: #otherwise, use the default name of "Unnamed Road"
            print("Travelling on Unnamed Road")

        # if the next road segment is on a bus route
        if (route[i], route[i+1]) in road_segments:
            # if the previous road segment is on a bus route
            if (route[i-1], route[i]) in road_segments:
                if set(road_segments[(route[i-1], route[i])]).intersection(set(road_segments[(route[i], route[i+1])])) != {}: # if possible, stay on the same bus
                    # you are on a bus, so it will drain battery
                    # Future Work: devise a more sophisticated battery draining/charging formula
                    battery -= (0.005 * x[0]["length"])
                    print("%.2f%% battery remaining" % battery)
                    print("Remain on bus " + prev_route)
                else: # the previous and next roads are both on bus routes, but have no common bus routes
                    print("Switch to bus " + curr_route)
                    battery = 100 # if getting on a new bus, assume it is initially at full charge
            else: # previous road was not on a bus route
                print("Get on bus " + curr_route)
                battery = 100 # if getting on a new bus, assume it is initially at full charge
                # Future Work: determine the battery charge of any bus atany given point
        else: # not on a bus route
            if (route[i-1], route[i]) in road_segments: #if the previous route is on a bus route but the next one is not
                print("Get off bus " + prev_route)

# A method to add the coordinates for each node in each route
def add_points(nodes):
    points = []
    for i in range(len(nodes)-1):
        a = nodes[i]
        b = nodes[i+1]
        length, path = nx.bidirectional_dijkstra(G, a, b, 'length')
        points += path # find the shortest path between each pair of nodes, and add these to the list of points
    points2 = points[:1] # remove duplicates, as these cause errors in plotting, as no node is connected to itself
    for j in points:
        if j != points2[-1]:
            points2.append(j)
    points3 = []
    for k in points2:
        points3.append([G.nodes[k]['y'], G.nodes[k]['x']]) # add the coordinates of each point

    # return the list of nodes and the list of coordinates
    return points2, points3


# colours every segment of road on each bus route in the colour(s) of its route(s)
def colour_lines(points, routenum, routemap):

    for i in range(len(points) - 1): #each segment is two points, so the last segment will start with the second-last stop
        if (tuple(points[i]), tuple(points[i+1])) not in road_segments:
            road_segments[(tuple(points[i]), tuple(points[i+1]))] = [routenum]
        else: #this segment is on another route, so add this route to the list
            if routenum not in road_segments[(tuple(points[i]), tuple(points[i+1]))]:
                road_segments[(tuple(points[i]), tuple(points[i+1]))] += [routenum]
            #otherwise, the segment is already on this route, don't re-add it


# add the label and place it at the bottom of the GUI
label = tk.Label(top, text="Select a start point by clicking on the map above")
label.pack(side=tk.BOTTOM)

# add instructions for how to zoom and pan to the GUI
instructions = tk.Label(top, text="Click the four-headed arrow to toggle panning and zooming, if it is on, click and drag in any direction with the left mouse to pan, or with the right mouse to zoom. Click outside the map when done to confirm choice.")
instructions.pack(side=tk.BOTTOM)

# convert the Excel file to a format that is readable by Python
file = pd.ExcelFile(timetable)

# create the shortest path and all bus routes, and add them to the map given the start and end points
# parameters: start = user's selected start point
# end = user's selected end point
# routemap = the map to place the routes on
def setup(start, end, routemap):

    # a list of all points on every route
    all_points = []

    n = 0 #unique ID for each bus stop, used to find nearest bus stop to start and end points

    # for each route
    for num in route_dict.keys():
        # set the route's colour
        linecolour = route_dict[num]["line_colour"]
        
        # read the bus stop data for that route into a data frame
        df = pd.read_excel(file, route_dict[num]["stops"], header=None)
        print("Excel file read")
        
        # get the stop coordinates from the data frame and add the coordinates to the list of stops
        stops = []
        stops_coords = df.iloc[:, 2]
        add_coords(stops_coords, stops)
        
        # add each bus stop to the list of all bus stops, and record the coordinates
        nodes = []
        xs = []
        ys = []
        for p in stops:
            all_bus_stops.append((n, p))
            xs += [p[1]]
            ys += [p[0]]
            n += 1
        
        # find the nearest edges for each edge on the route
        nearest = ox.get_nearest_edges(G, xs, ys, 'balltree')

        # 
        for s in range(len(nearest)):
            nodes += [G.nodes()[nearest[s][0]]['osmid']]
            if s == len(nearest) - 1:
                nodes += [G.nodes()[nearest[s][1]]['osmid']]
        
        # add all bus stops and charging stations to the map
        add_markers(stops, num, charging_stations)
        
        # add all points on each route
        points, points_2 = add_points(nodes)
        
        # colour in each bus route
        colour_lines(points_2, num, routemap)
        print("Points added")
        all_points.append(points)
        
        # add each feature group to the map to allow toggling
        routemap.add_child(fg_dict[num])

    for j in road_segments.keys(): # for each road segment
        for k in range(len(road_segments[j])): # for each route that this segment is part of
            current = road_segments[j][k]
            dashes = "%s %s" % (10 * (len(road_segments[j]) - k), 10 * k) # splits the line into n equal parts, no matter how big n is
            pl = f.PolyLine(j, color=route_dict[current]["line_colour"], dash_array=dashes, weight=10)
            fg_dict[current].add_child(pl)

    # find the nearest nodes to the start and end selected by the user
    start_node = ox.get_nearest_node(G, start)
    end_node = ox.get_nearest_node(G, end)

    # add markers to the user's chosen start and end
    f.Marker(location=start, icon=f.Icon(color="green")).add_to(routemap)
    f.Marker(location=end, icon=f.Icon(color="green")).add_to(routemap)
    print("Start and end point markers added")

    # find nearest bus stop to start and end points, use these for route in display_route, then use double Dijkstra for walking from start to first stop, and from last stop to end
    start_stop = get_nearest_bus_stop(G, G.nodes()[start_node]['y'], G.nodes()[start_node]['x'])
    end_stop = get_nearest_bus_stop(G, G.nodes()[end_node]['y'], G.nodes()[end_node]['x'])

    # find the first and last bus stops on the route
    start_stop_node = ox.get_nearest_node(G, all_bus_stops[start_stop][1])
    end_stop_node = ox.get_nearest_node(G, all_bus_stops[end_stop][1])

    # a list of all road segments that are on bus routes
    bus_routes = []

    for edge in G.edges():
        # for each edge in the graph, check if it is on a bus route
        start = (G.nodes()[edge[0]]['y'], G.nodes()[edge[0]]['x'])
        end = (G.nodes()[edge[1]]['y'], G.nodes()[edge[1]]['x'])
        edge_road = G[edge[0]][edge[1]]
        if (start, end) in road_segments: #bus routes should be prioritized, as bussing takes less time than walking
            bus_routes += [(edge[0], edge[1])]
        edge_road[0]['time'] = edge_road[0]['length'] / WALKING_SPEED # all non-bus routes take longer to traverse as they must be walked

    for edge in bus_routes: # all bus routes get a smaller time value because they take less time to traverse
        # the edges are one-way as the bus cannot travel in reverse, and bus routes are not the same going forward as they are going back
        G.add_edge(edge[0], edge[1], oneway=True, length=edge_road[0]['length'], time=edge_road[0]['length']/BUS_SPEED)

    #compute path from start to first stop, and last stop to end, and add these to path
    before_path = nx.bidirectional_dijkstra(G, start_node, start_stop_node, 'time')[1] #walking up to first stop
    temp_path = nx.bidirectional_dijkstra(G, start_stop_node, end_stop_node, 'time')[1]
    after_path = nx.bidirectional_dijkstra(G, end_stop_node, end_node, 'time')[1] #walking from last stop
    
    # computing the three parts of the path
    xs = []
    ys = []
    path = []
    for a in before_path:
        b = G.nodes()[a]
        path += [(b['y'], b['x'])]
    for p in temp_path:
        q = G.nodes()[p]
        path += [(q['y'], q['x'])]
    for s in after_path:
        t = G.nodes()[s]
        path += [(t['y'], t['x'])]

    # add the shortest path to the map and the legend
    fg = f.FeatureGroup(name = lgd_txt.format(txt="Shortest Path", col="black"))
    pl = f.PolyLine(path, color="black", weight=4) #shortest path goes over bus routes, must be narrower
    fg.add_child(pl)
    routemap.add_child(fg)

    # display the shortest path along with directions to the user
    display_route(path)

    # add the legend to the map, and save the map as a HTML file
    f.map.LayerControl('topleft', collapsed= False).add_to(routemap)
    routemap.save('route.html')
    print("Route saved")
    exit()

# get the OSM tiles for the initial selection map
tiles.OSM.get_image = image_spoof
osm_img = tiles.OSM()

# plot the map on a graph using Matplotlib
fig = plt.figure(figsize=(12,6))
ax = plt.axes(projection=ccrs.PlateCarree()) #creates a map with standard lat-long coordinates

#to use a map of another area, put in min longitude, max longitude, min latitude and max latitude
extent = [-8.65, -8.33, 51.7, 51.96] #area of Cork covered by bus routes
ax.set_extent(extent)

scale = 15 # higher scales consume more memory, especially for large areas, so 15 is the maximum value that will work
ax.add_image(osm_img, int(scale))

# used to store the coordinates of the user's clicks on the map
coords = []

MAX_CLICK_LENGTH = 0.1 #used to distinguish click-and-drag for zooming/panning from clicking to choose a point

# functions for handling user clicks
def on_click(event, ax):
    ax.time_onclick = time.time()

def on_release(event, ax):
    if event.inaxes == ax: #if click is within the axes
        if event.button == 1 and ((time.time() - ax.time_onclick) < MAX_CLICK_LENGTH): #only count a single left click, not click and drag or right click
            coords.append((event.ydata, event.xdata))
            print(event.xdata, event.ydata)
            if len(coords) == 1:
                label['text'] = "Select an end point by clicking on the map above"
            if len(coords) == 2:
                canvas.mpl_disconnect(c1)
                canvas.mpl_disconnect(c2)
                setup(coords[0], coords[1], map_of_cork)

# add in canvas to display map
canvas = FigureCanvasTkAgg(fig, master=top)
canvas.get_tk_widget().pack(side=tk.TOP)

# add event handlers to click-based events
c1 = canvas.mpl_connect('button_press_event', lambda event: on_click(event, ax))
c2 = canvas.mpl_connect('button_release_event', lambda event: on_release(event, ax))

# add toolbar to GUI
toolbar = NavigationToolbar2Tk(canvas, top)
canvas._tkcanvas.pack(side=tk.TOP)

# add quit button to GUI
button = tk.Button(master=top, text="Quit", command=exit)
button.pack(side=tk.BOTTOM)

top.mainloop()
