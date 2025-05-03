# Routing with RSMT

## Overview

Routing is a critical step in the physical synthesis process that involves connecting pins with metal wires while minimizing total net length. In this repo, I use **Sequential Steiner Tree Heuristic Algorithm** to construct the Rectilinear Steiner Minimum tree and connect all pins. By sequentially linking each pin based on the shortest-first principle, this algorithm efficiently reduces the net length while assuring the connectivity of all pins. Also, I remove the redundant edges from net loops caused by net overlap by depth-first searching (DFS).

My implementation uses Python 3.10.16 and only imports minimal libraries like pyplot and math. No hyperparameter needs to be specified in this repo. I tested the algorithm on four benchmark problems of increasing size: 10, 50, 150, and 200 pins. The final results are listed as follows (`verify()` function has ensured **all pins are connected**):

| #Pins | Wirelength | Runtime (s) |
|-------|------------|-------------|
| 10    | 308.0      | 0.20        |
| 50    | 3531.0     | 0.24        |
| 150   | 16086.0    | 0.56        |
| 200   | 25981.0    | 0.91        |

## Usage

After changing to `rsmt_generation` directory, simply run the command `python rsmt_generation.py` to execute the program. The program will generate output files and visualization figs in `./output` and `./fig`. The input traces have already been provided in the directory `./3-RSMT-v0`.