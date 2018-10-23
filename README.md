# 20181008_Bioinformatics-Tools
## Find superstring from sequencing reads
### Algorithm
- This is a **greedy** Hamiltonian Path problem. More info see here:
- http://mathworld.wolfram.com/HamiltonianPath.html
### What does my code do
- The file 'Find_supertring.py' allows you to assemble the sequencing reads to superstring.
- It first add the reads to a graph, then based on the overlaps between reads, the edges are added.
- The superstring are traversed along the edge.
## Affine gap global alignment with context-sensitive gap costs
### Algorithm
```
M(i, j) = max{ M(i - 1, j - 1) + S(x_i, y_j), I_x(i - 1, j - 1) + S(x_i, y_j), I_y(i - 1, j - 1) + S(x_i, y_j) }
I_x(i, j) = max{ M(i - 1, j) + g(y_j, y_{j+1}) + s, I_x(i - 1, j) + s }
I_y(i, j) = max{ M(i, j - 1) + g(x_i, x_{i+1}) + s, I_y(i, j - 1) + s }
```
### What does my code do
- Perform globle sequencing algiment between two sequences.
- Use affign gap penalty *(gap opening penalty and gap extension penalty)* 
- The gap costs are dependent on the base pairs before and after the gap
