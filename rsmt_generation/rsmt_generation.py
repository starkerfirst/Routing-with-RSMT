################################################
#ELEC6910D: RSMT generation
#Author: YANG Bohan
#Algorithm: Sequential Steiner Tree Heuristic
################################################

import os
import matplotlib.pyplot as plt
import math
import time

verbose = False
input_dir = '3-RSMT-v0'

class Edge:
    def __init__(self, p1, p2):
        assert(p1[0] == p2[0] or p1[1] == p2[1])
        assert(isinstance(p1, tuple) and isinstance(p2, tuple))
        self.direction = "vertical" if p1[0] == p2[0] else "horizontal"
        self.invalid = False
        self.visited = False

        if p1 == p2:
            self.invalid = True

        if self.direction == "vertical":
            if p1[1] < p2[1]:
                self.p1 = p1
                self.p2 = p2
            else:
                self.p1 = p2
                self.p2 = p1
        else:
            if p1[0] < p2[0]:
                self.p1 = p1
                self.p2 = p2
            else:
                self.p1 = p2
                self.p2 = p1

    def string(self):
        return f"{self.p1[0]} {self.p1[1]} {self.p2[0]} {self.p2[1]}"
    
    def get_all_points(self):
        if self.direction == "vertical":
            return [(self.p1[0], y) for y in range(min(self.p1[1], self.p2[1]), max(self.p1[1], self.p2[1]) + 1)]
        else:
            return [(x, self.p1[1]) for x in range(min(self.p1[0], self.p2[0]), max(self.p1[0], self.p2[0]) + 1)]
        
    def point_is_on_edge(self, point):
        if self.direction == "vertical":
            return point[0] == self.p1[0] and min(self.p1[1], self.p2[1]) <= point[1] <= max(self.p1[1], self.p2[1])
        else:
            return point[1] == self.p1[1] and min(self.p1[0], self.p2[0]) <= point[0] <= max(self.p1[0], self.p2[0])
        
    # Check if this edge is contained in another edge.
    def is_contained_by(self, other):
        if self.direction == "vertical" and other.direction == "vertical":
            return (self.p1[0] == other.p1[0] and 
                    self.p1[1] >= other.p1[1] and 
                    self.p2[1] <= other.p2[1])
        elif self.direction == "horizontal" and other.direction == "horizontal":
            return (self.p1[1] == other.p1[1] and 
                    self.p1[0] >= other.p1[0] and 
                    self.p2[0] <= other.p2[0])
        return False
    
    def has_34port_overlap(self, other):
        if self.p1 == other.p1 or self.p1 == other.p2 or self.p2 == other.p1 or self.p2 == other.p2:
            return False
        if self.direction == "vertical" and other.direction == "horizontal":
            # Check if the horizontal edge crosses the vertical edge
            return (other.p1[0] <= self.p1[0] <= other.p2[0] and 
                    self.p1[1] <= other.p1[1] <= self.p2[1])
        elif self.direction == "horizontal" and other.direction == "vertical":
            # Check if the vertical edge crosses the horizontal edge
            return (other.p1[1] <= self.p1[1] <= other.p2[1] and
                    self.p1[0] <= other.p1[0] <= self.p2[0])
        return False

    def length(self):
        assert(not self.invalid)
        return math.dist(self.p1, self.p2)
    

class RSMT:
    def __init__(self, numPins, loc_list):
        self.num_pins = numPins
        self.pin_loc_list = loc_list
        self.start_point = None
        assert(self.num_pins == len(self.pin_loc_list))

        # Initialize other necessary attributes
        self.tree_edges = []
        self.pruned_tree_edges = []
        self.steiner_points = []
        self.remaining_nodes = list(range(self.num_pins))
        
    def get_wirelength(self):
        total_wirelength = 0
        for edge in self.tree_edges:
            total_wirelength += edge.length()
        return total_wirelength
    
    def print_tree(self):
        for edge in self.tree_edges:
            print(edge.string())

    def clear_edge_state(self):
        for edge in self.tree_edges:
            edge.visited = False

    def remove_edge(self, p1, p2):
        remove_list = []
        for edge in self.tree_edges:
            if edge.p1 == p1 and edge.p2 == p2 and (p1, p2) not in remove_list:
                remove_list.append((p1, p2))
            elif edge.p1 == p2 and edge.p2 == p1 and (p2, p1) not in remove_list:
                remove_list.append((p2, p1))
        for edge in self.tree_edges:
            if (edge.p1, edge.p2) in remove_list or (edge.p2, edge.p1) in remove_list:
                self.tree_edges.remove(edge)

    def find_neighbors(self, point, except_neighbor=None):
        # find the neighbors of point
        neighbors = []
        for edge in self.tree_edges:
            assert(not edge.invalid)
            if edge.p1 == point and edge.p2 != except_neighbor:
                neighbors.append(edge.p2)
            elif edge.p2 == point and edge.p1 != except_neighbor:
                neighbors.append(edge.p1)
        return neighbors
    
    def remove_invalid_edges(self):
        # remove the invalid edges
        tree_edges = []
        for edge in self.tree_edges:
            if not edge.invalid:
                tree_edges.append(edge)
        self.tree_edges = tree_edges

    # Sequential Steiner Tree Heuristic algorithm
    def compute_rsmt(self):
        assert(self.tree_edges == [])
        assert(self.steiner_points == [])

        # create the first point pair
        min_dist = float('inf')
        for i in range(self.num_pins):
            for j in range(i+1, self.num_pins):
                dist = math.dist(self.pin_loc_list[i], self.pin_loc_list[j])
                if dist < min_dist:
                    min_dist = dist
                    p1, p2 = i, j
        self.remaining_nodes.remove(p1)
        self.remaining_nodes.remove(p2)
        self.start_point = self.pin_loc_list[p1]
        # breakpoint()

        # construct the minimum bounding box
        L_shape_1 = (Edge(self.pin_loc_list[p1], (self.pin_loc_list[p2][0], self.pin_loc_list[p1][1])), 
                            Edge((self.pin_loc_list[p2][0], self.pin_loc_list[p1][1]), self.pin_loc_list[p2]))
        L_shape_2 = (Edge(self.pin_loc_list[p1], (self.pin_loc_list[p1][0], self.pin_loc_list[p2][1])),
                            Edge((self.pin_loc_list[p1][0], self.pin_loc_list[p2][1]), self.pin_loc_list[p2]))
        minimum_bounding_box = (L_shape_1, L_shape_2)
        p_mbb_list = []
        for L_shape in minimum_bounding_box:
            for edge in L_shape:
                for point in edge.get_all_points():
                    p_mbb_list.append(point)

        while len(self.remaining_nodes) > 0:
            # print(f"Remaining nodes: {self.remaining_nodes}")
            # breakpoint()
                        
            # find the closest point pair (p_mbb_list, p_c)
            min_dist = float('inf')
            for i in range(len(self.remaining_nodes)):
                for j in range(len(p_mbb_list)):
                    dist = math.dist(self.pin_loc_list[self.remaining_nodes[i]], p_mbb_list[j])
                    if dist < min_dist:
                        min_dist = dist
                        p1, p2 = i, j # p2 is p_mbb
            p1 = self.remaining_nodes.pop(p1)
            

            # remove the L_shape which p is not on
            is_found = False
            for L_shape in minimum_bounding_box:
                for edge in L_shape:
                    if edge.point_is_on_edge(p_mbb_list[p2]):
                        is_found = True
                        self.tree_edges.append(L_shape[0])
                        self.tree_edges.append(L_shape[1])
                        break
                # if point is on both edges, we choose the first one
                if is_found:
                    break
            assert(is_found)

            # construct the minimum bounding box
            L_shape_1 = (Edge(self.pin_loc_list[p1], (p_mbb_list[p2][0], self.pin_loc_list[p1][1])), 
                                Edge((p_mbb_list[p2][0], self.pin_loc_list[p1][1]), p_mbb_list[p2]))
            L_shape_2 = (Edge(self.pin_loc_list[p1], (self.pin_loc_list[p1][0], p_mbb_list[p2][1])),
                                Edge((self.pin_loc_list[p1][0], p_mbb_list[p2][1]), p_mbb_list[p2]))
            minimum_bounding_box = (L_shape_1, L_shape_2)
            p_mbb_list = []
            for L_shape in minimum_bounding_box:
                for edge in L_shape:
                    for point in edge.get_all_points():
                        p_mbb_list.append(point)
        
        # add the last edge (first L_shape)
        self.tree_edges.append(minimum_bounding_box[0][0])
        self.tree_edges.append(minimum_bounding_box[0][1])

        # prune the redundant edges
        self.edge_redundancy_process()
        
    # prune the redundant/overlap edges
    # generate steiner points
    def edge_redundancy_process(self):
        self.clear_edge_state()

        # break long edges with overlap into short edges
        is_converged = False
        while not is_converged:
            is_converged = True
            for i in range(len(self.tree_edges)):
                for j in range(i+1, len(self.tree_edges)):
                    edge_i = self.tree_edges[i]
                    edge_j = self.tree_edges[j]
                    if edge_i.has_34port_overlap(edge_j):
                        is_converged = False
                        # find the intersection point
                        if edge_i.direction == "vertical":
                            junction = (edge_i.p1[0], edge_j.p1[1])
                        else:
                            junction = (edge_j.p1[0], edge_i.p1[1])
                        new_edge_1 = Edge(edge_i.p1, junction)
                        new_edge_2 = Edge(junction, edge_i.p2)
                        new_edge_3 = Edge(edge_j.p2, junction)
                        new_edge_4 = Edge(junction, edge_j.p1)
                        self.tree_edges.remove(edge_i)
                        self.tree_edges.remove(edge_j)
                        self.tree_edges.append(new_edge_1)
                        self.tree_edges.append(new_edge_2)
                        self.tree_edges.append(new_edge_3)
                        self.tree_edges.append(new_edge_4)
                    
        # remove the invalid edges
        self.remove_invalid_edges()
        
        # Start DFS from the start point
        self.visited = []
        self.dfs(self.start_point, None)
        self.verify()

        # remove steiner points with no more than 2 neighbors
        for point in self.steiner_points:
            neighbors = self.find_neighbors(point)
            if len(neighbors) <= 2:
                self.steiner_points.remove(point)

        print(f"Steiner points: {self.steiner_points}")
        
    # Recursive DFS to find Steiner points and remove redundant edges.
    def dfs(self, current, parent):
        # print(f"Visiting: {current}")
        if current in self.visited:
            return
        self.visited.append(current)

        # if current == (228, 36):
        #     breakpoint()
        
        neighbors = self.find_neighbors(current, except_neighbor=parent)
        remove_list = []
        for neighbor in neighbors:
            if neighbor in self.visited:
                # there is a cycle
                remove_list.append(neighbor)
            
        for neighbor in remove_list:    
            self.remove_edge(current, neighbor)
            neighbors.remove(neighbor)

        # review the neighbors
        neighbors = self.find_neighbors(current)
        if len(neighbors) == 1 and current not in self.pin_loc_list:
            self.remove_edge(parent, current)
            return

        if current not in self.pin_loc_list and current not in self.steiner_points and len(neighbors) > 2:
            self.steiner_points.append(current)

        # Recursively visit all neighbors
        neighbors = self.find_neighbors(current, except_neighbor=parent)
        neighbors.sort(key=lambda x: math.dist(current, x))
        for neighbor in neighbors:
            self.dfs(neighbor, current)
        
    # connectivity check
    def verify(self):
        for pin in self.pin_loc_list:
            if pin not in self.visited:
                print(f"Pin {pin} is not connected to the tree.")
                raise ValueError("Tree is not connected.")

def collect_data(file):
    loc_list = []
    file_path = os.path.join(input_dir, file)
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Parse first line: Net <net_id> <num_pins>
        first_line = lines[0].strip().split()
        if len(first_line) >= 3 and first_line[0] == "Net":
            net_id = int(first_line[1])
            num_pins = int(first_line[2])
        else:
            raise ValueError(f"Invalid input format in file {file}")
        
        # Read each pin's coordinates
        x_min = float('inf')
        x_max = 0
        y_min = float('inf')
        y_max = 0
        for i in range(1, num_pins + 1):
            pin_data = lines[i].strip().split()
            if len(pin_data) >= 3:
                x_coord = int(pin_data[1])
                y_coord = int(pin_data[2])
                loc_list.append((x_coord, y_coord))
                if x_coord < x_min:
                    x_min = x_coord
                if x_coord > x_max:
                    x_max = x_coord
                if y_coord < y_min:
                    y_min = y_coord
                if y_coord > y_max:
                    y_max = y_coord
                
        x_range = (x_min, x_max)
        y_range = (y_min, y_max)
        return num_pins, loc_list, x_range, y_range

def plot_solution(rsmt, file, x_range, y_range):
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot edges first (bottom layer)
    for edge in rsmt.tree_edges:
        p1, p2 = edge.p1, edge.p2
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='green', linewidth=1, zorder=1)
    
    # Plot pins (middle layer)
    x_coords = [loc[0] for loc in rsmt.pin_loc_list]
    y_coords = [loc[1] for loc in rsmt.pin_loc_list]
    ax.scatter(x_coords, y_coords, color='blue', s=5, label='Pins', zorder=2)
    
    # Plot pin labels
    # for i, (x, y) in enumerate(rsmt.pin_loc_list):
    #     ax.annotate(str(i), (x, y), xytext=(5, 5), textcoords='offset points')
    
    # Plot Steiner points (top layer)
    if rsmt.steiner_points:
        steiner_x = [p[0] for p in rsmt.steiner_points]
        steiner_y = [p[1] for p in rsmt.steiner_points]
        ax.scatter(steiner_x, steiner_y, color='red', s=2, label='Steiner Points', zorder=3)
    
    ax.set_title(f'RSMT Solution')
    ax.legend()
    # ax.grid(True)

    if not os.path.exists('fig'):
        os.makedirs('fig')
    
    plt.savefig(f'fig/rsmt_{file}.png', dpi=300)

def export_solution(rsmt, file_name):
    num_folder = file_name.split('.')[0]
    
    output_folder = f"output/{num_folder}"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    steiner_file_path = os.path.join(output_folder, f'steiner_points_list.txt')
    with open(steiner_file_path, 'w') as f:
        for point in rsmt.steiner_points:
            f.write(f"{point[0]} {point[1]}\n")
    
    edge_file_path = os.path.join(output_folder, f'edge_list.txt')
    with open(edge_file_path, 'w') as f:
        for edge in rsmt.tree_edges:
            f.write(f"{edge.p1[0]} {edge.p1[1]} {edge.p2[0]} {edge.p2[1]}\n")
    
    print(f"RSMT solution exported to {output_folder}")

def main():
    file_list = os.listdir(input_dir)
    file_list = sorted([file for file in file_list if file.endswith('.txt')], key=lambda x: int(x.split('.')[0]))

    for file in file_list:
        start = time.time()
        print(f'Using file: {file}')
        num_pins, loc_list, x_range, y_range = collect_data(file)
        print(f"Number of pins: {num_pins}")
        print(f"X range: {x_range}")
        print(f"Y range: {y_range}")

        # Initialize RSMT
        rsmt = RSMT(num_pins, loc_list)
        rsmt.compute_rsmt()
        rsmt.verify()

        # Print results
        print(f"Total wirelength: {rsmt.get_wirelength()}")

        export_solution(rsmt, file)
        plot_solution(rsmt, file.split('.')[0], x_range, y_range)

        end = time.time()
        print(f'Time taken: {end - start:.2f} seconds')
        # breakpoint()

if __name__ == "__main__":
    main()