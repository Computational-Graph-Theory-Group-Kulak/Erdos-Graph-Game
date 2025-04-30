# Erdos-Graph-Game
This repository contains a solver for a graph-colouring game designed by Erdös which is created for the article "On edge-colouring-games by Erdös, and Bensail and Mc Inerney".

The latest version of this program can be obtained from <https://github.com/Computational-Graph-Theory-Group-Kulak/Erdos-Graph-Game>.

This program can be used to find the result of an optimal game on a clique or Colex-graph of a specified order. It makes use of datastructures and methods from [`nauty`](https://pallini.di.uniroma1.it/) and [`K2-Hamiltonian Graphs`](https://github.com/JarneRenders/K2-Hamiltonian-Graphs).

### Installation

This requires a working shell and `make`.

- Download, extract and configure [`nauty`](https://pallini.di.uniroma1.it/).
- Configure and compile the nauty libraries for multithreading using: `./configure --enable-tls` and `make`.
- Compile using: 
	* `make 64-bit` to create a binary for the 64 bit version
    * `make 128-bit` to create a binary for the 64 bit version

The 64 bit version is significantly faster than the 128 bit version, and it is advised to use this version 


### Usage of erdos-solver

This helptext can be found by executing `./erdos-solver -h`.

Usage: `./erdos-solver n [-h] [-g] [-t] [-b] [-p] red-graph blue-graph `

When the -g flag is not set, graphs are read from stdin in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>. When the -g flag is set graphs are read as bitsets from stdin.

The order of the graphs is always the first argument. The blue startgraph argument should always appear after the red startgraph argument. Otherwise, the order in which the arguments appear does not matter. Be careful not to put an argument immediately after one with an option. E.g. -g#b will not recognise the -b argument.

Mandatory arguments to long options are mandatory for short options too.
```
  -b, --bias=#              blue can select # more edges than red in its turn
  -t, --threads=#           use a maximum of # threads 
  -g, --generator-used      the input graphs have been provided by the generator in bitset format
  -p, --starting-player=#   the starting player in blue if #==2 or red otherwise
  -h, --help                print help message
  n                         the order of the graph on which the game is played
  red-graph                 the start graph of the red player in graph6 format
  blue-graph                the start graph of the blue player in graph6 format
```

### Usage of gen-colex

This helptext can be found by executing `./gen-colex -h`.

Usage: `./gen-colex n [-h] [-b] [-p] base-graph red-graph blue-graph `

When the -g flag is not set, graphs are read from stdin in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>. When the -g flag is set graphs are read as bitsets from stdin.

The order of the graphs is always the first argument. The blue startgraph argument should always appear after the red startgraph argument. both startgraphs should always appear after the base graph. Otherwise, the order in which the arguments appear does not matter. Be careful not to put an argument immediately after one with an option. E.g. -p#b will not recognise the -b argument.

Mandatory arguments to long options are mandatory for short options too.
```
  -b, --bias=#              blue can select # more edges than red in its turn
  -p, --starting-player=#   the starting player in blue if #==2 or red otherwise
  -h, --help                print help message
  n                         the order of the graph on which the game is played
  base-graph                the base graph with all allowed edges in graph6 format
  red-graph                 the start graph of the red player in graph6 format
  blue-graph                the start graph of the blue player in graph6 format
```

### Examples
`bash Erdos-Game-Generic.sh 7 F~~~w F???? F???? 0 1 1`
Plays the game using our generator on K_7 with an empty start graph, 0 bias between the two players, 1 thread for the solver and red starts (in that order)

`bash Erdos-Game-Cliques.sh 7`
Plays the game using geng on K_7 with an empty start graph, 0 bias between the two players, 1 thread for the solver and red starts
