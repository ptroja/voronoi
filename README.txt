This utility computes a safe the path for a mobile robot on a map.

It computes a Voronoi diagram [1] of a map using Qhull [2].
Then vertices at start and goal positions are added to the diagram.
Shortest path is computed with Dijkstra's algorithm [3] using Boost Graph
Library implementation [4]. Finally path is simplified with a Douglas–Peucker
algorithm [5], based on [6] implementation.

Map images and display is handled with a CImg library [7] (tested with both
CImg-1.3.9 and CImg-1.4.4).

Left-mouse button marks a start position on a map, while right-mouse button
marks a goal position.

USAGE:
./voronoi map_image_file [minFacetArea] [additional_Qhull_arguments] ...

TODO:
1. Voronoi diagram as computed with Qhull is not perfect/complete. Voronoi
vertices are detected between map pixels and filtering of this case is not
perfect. The density of a diagram can be tuned with minFacetArea parameter.

You can observe this sub-pixel vertices with the following Matlab commands:
> cave = imread('cave.png');
> invcave = (cave == 0);
> [x,y] = find (invcave);
> voronoi(x,y);

Probably this can be solved using Fortune's algorithm for computing Voronoi
diagram [7] (link to the C-implementation included). 

2. The minimum distance to the obstacles is not verified. This can be done
at the beginning by enlarging occupied areas on the map.

3. After #ifdef'ing CImg routines this should be implemented into
a Player/Stage [9] navigation driver.

AUTHOR:
Piotr Trojanek <piotr.trojanek (at) gmail.com>  

[1] http://en.wikipedia.org/wiki/Voronoi_diagram
[2] http://www.qhull.org/
[3] http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
[4] http://www.boost.org/doc/libs/
[5] http://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm
[6] http://softsurfer.com/Archive/algorithm_0205/algorithm_0205.htm
[7] http://cimg.sourceforge.net/
[8] http://en.wikipedia.org/wiki/Fortune%27s_algorithm
[9] http://playerstage.sourceforge.net/