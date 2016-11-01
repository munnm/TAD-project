# This is based on depth first search

# An outer loop scans all nodes of graph and starts a search from every node

# Node neighbors are added to the cycle path

# Recursion ends if no more non-visited neighbors can be added

# A new cycle is found if the path is longer than 4 nodes and next neighbor is the start of the path

# To avoid duplicate cycles, the cycles are normalized by rotating the smallest node to the start

import time

#


graph = [[1,2],[1,3],[1,4],[2,3],[3,4],[2,6],[4,6],[7,8],[8,9],[9,7]]
graph = [[1, 2], [3, 4], [3, 5], [6, 10], [7, 9], [7, 10], [8, 9], [12, 13], [15, 16], [17, 18], [19, 20], [21, 22], [23, 24], [23, 25], [27, 28], [30, 31], [30, 32], [34, 35], [36, 37], [38, 39], [44, 46], [44, 47], [45, 46], [48, 49], [48, 50], [51, 52], [51, 53], [54, 55], [56, 57], [56, 58], [60, 61], [63, 64], [69, 70], [71, 72], [76, 77], [78, 79], [83, 85], [84, 85], [86, 87], [89, 91], [89, 92], [90, 91], [93, 94], [93, 96], [95, 96], [98, 99], [100, 101], [102, 103], [105, 106], [107, 108], [112, 113], [116, 117], [121, 122], [125, 126], [128, 129], [128, 130], [131, 132], [136, 137], [139, 140], [139, 141], [146, 147], [149, 150], [151, 152], [158, 159], [158, 160], [161, 162], [163, 164], [165, 166], [167, 168], [170, 171], [172, 173], [178, 179], [178, 180], [182, 183], [184, 185], [187, 188], [190, 192], [191, 192], [194, 195], [196, 197], [200, 201], [200, 202], [206, 207], [208, 209], [213, 214], [213, 215], [217, 218], [219, 221], [220, 221], [222, 223], [225, 226], [227, 228], [229, 230], [229, 231], [234, 235], [240, 241], [242, 243], [242, 244], [248, 249], [250, 251], [250, 252], [253, 255], [254, 255], [256, 258], [257, 258], [259, 262], [260, 262], [261, 262], [263, 264], [265, 266], [267, 268], [269, 270], [272, 273], [274, 275], [277, 279], [278, 279], [281, 282], [284, 285], [286, 288], [287, 288], [289, 290], [289, 292], [291, 292], [293, 296], [294, 295], [294, 296], [297, 298], [299, 300], [304, 305], [307, 308], [309, 310], [311, 312], [313, 314], [313, 315], [316, 317], [318, 319], [321, 322], [324, 327], [325, 327], [326, 327], [328, 329], [335, 336], [337, 338], [339, 340], [341, 342], [341, 343], [348, 350], [349, 350], [352, 353], [354, 355], [356, 358], [357, 358], [359, 360], [366, 367], [368, 369], [370, 371], [372, 373], [374, 375], [377, 378], [387, 388], [389, 390], [389, 391], [393, 394]]
cycles = []


def main():
    global graph
    global cycles
    start_time = time.time()
    for edge in graph:
      for vertex in edge:
        findNewCycles([vertex])
    print(cycles)
    print("---%s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    for edge in graph:
      for vertex in edge:
        findNewCycles_revised([vertex])
    print(cycles)
    print("---%s seconds ---" % (time.time() - start_time))
    
    
    
def findNewCycles_revised(path):
    start_node = path[0]
    next_node = None
    sub = []

    for edge in graph:
        node1, node2 = edge
        # stop search once you're too far along in the graph
        if node1 <= path[0]:
            if start_node in edge:
                    if node1 == start_node:
                        next_node = node2
                    else:
                        next_node = node1
            if not visited(next_node, path):
		# neighbor node not on path yet
                    sub = [next_node]
                    sub.extend(path)
		# explore extended path
                    findNewCycles(sub);
            elif len(path) > 3  and next_node == path[-1]:
		# cycle found
                    p = rotate_to_smallest(path);
                    inv = invert(p)
                    if isNew(p) and isNew(inv):
                        cycles.append(p)

def findNewCycles(path):
    start_node = path[0]
    next_node = None
    sub = []

    for edge in graph:
        node1, node2 = edge
        if start_node in edge:
                if node1 == start_node:
                    next_node = node2
                else:
                    next_node = node1
        if not visited(next_node, path):
		# neighbor node not on path yet
                sub = [next_node]
                sub.extend(path)
		# explore extended path
                findNewCycles(sub);
        elif len(path) > 3  and next_node == path[-1]:
		# cycle found
                p = rotate_to_smallest(path);
                inv = invert(p)
                if isNew(p) and isNew(inv):
                    cycles.append(p)


# next two function is for formating purpose only
def invert(path):
    return rotate_to_smallest(path[::-1])



# rotate cycle path such that it begins with the smallest node
def rotate_to_smallest(path):
    n = path.index(min(path))
    return path[n:]+path[:n]



# check if it is a unique path
def isNew(path):
    return (not path in cycles)



# check if a vertex has been visited
def visited(vertex, path):
    return (vertex in path)

main()
