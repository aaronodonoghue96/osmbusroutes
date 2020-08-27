# osmbusroutes
This repository contains a program to find the shortest path between two points using bus routes where possible, using OpenStreetMap data for map generation, called bus_routes_shortest_path. 
The other file, timetable_cork.xlsx is a Microsoft Excel file detailing the bus routes of Cork, including the coordinates and names of each bus stop along each route.
These files must be in the same directory when downloaded in order for the program to find the bus route data necessary to function.

Instructions to run this code:
1. Download the latest version of Anaconda from anaconda.com, and set it up on your machine
2. Download OSMNX and cartopy using conda-forge (these libraries do not come with Anaconda by default)
3. Set up an environment for the code to run in using the Anaconda Prompt, using the following command: 
4. Activate your environment using the command: conda activate <name of environment>
5. Navigate to the directory where the files from this repository are installed, using the Anaconda Prompt, and type: python bus_routes_shortest_path.py.
  
This code can also be run using a Jupyter Notebook (also part of Anaconda). To do this, type "jupyter notebook" into the Anaconda Prompt, which will open Jupyter and display your current directory. Navigate to the directory with bus_routes_shortest_path.py, and create a new Jupyter Notebook. Add in the code from bus_routes_shortest_path.py to the notebook and press the play button to run the program via the notebook.
