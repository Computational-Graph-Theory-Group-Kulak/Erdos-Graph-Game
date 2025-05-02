# Erdos-Graph-Game
This repository contains a solver for a graph-colouring game designed by Erdös which is created for the article "On edge-colouring-games by Erdös, and Bensail and Mc Inerney".

The latest version of this program can be obtained from <https://github.com/Computational-Graph-Theory-Group-Kulak/Erdos-Graph-Game>.

This program can be used to find the result of an optimal game on a clique or Colex-graph of a specified order. It makes use of datastructures and methods from [`nauty`](https://pallini.di.uniroma1.it/) and [`K2-Hamiltonian Graphs`](https://github.com/JarneRenders/K2-Hamiltonian-Graphs).

### Installation

This requires a working shell and `make`.

- Download, extract and configure [`nauty`](https://pallini.di.uniroma1.it/).
- Configure and compile the nauty libraries for multithreading using: `./configure --enable-tls` and `make`.
- Compile using: 
	* `make all-64bit` to create a binary for the 64 bit version

### Usage of Erdos-Game-Generic

Usage: `bash Erdos-Game-Generic.sh n start-graph red-graph blue-graph threadnumber starting-player`

The order of the arguments here is important.
```
  n                         the order of the graph on which the game is played
  red-graph                 the start graph of the red player in graph6 format
  blue-graph                the start graph of the blue player in graph6 format
  threadnumber              use a maximum of threadnumber threads 
  starting-player           the starting player in blue if starting-player==2 or red otherwise
```

### Usage of Erdos-Game-Cliques

Usage: `bash Erdos-Game-Cliques.sh n`

```
  n                         the order of the graph on which the game is played
```


### Usage of erdos-solver

This helptext can be found by executing `./erdos-solver -h`.

Usage: `./erdos-solver n [-h] [-g] [-t] [-b] [-p] red-graph blue-graph `

When the -g flag is not set, graphs are read from stdin in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>. When the -g flag is set graphs are read as bitsets from stdin.

The order of the graphs is always the first argument. The blue startgraph argument should always appear after the red startgraph argument. Otherwise, the order in which the arguments appear does not matter. Be careful not to put an argument immediately after one with an option. E.g. -g#b will not recognise the -b argument.

```
  -b#                       blue can select # more edges than red in its turn
  -t#                       use a maximum of # threads 
  -g                        the input graphs have been provided by the generator in bitset format
  -p#                       the starting player in blue if #==2 or red otherwise
  -h                        print help message
  n                         the order of the graph on which the game is played
  red-graph                 the start graph of the red player in graph6 format
  blue-graph                the start graph of the blue player in graph6 format
```

### Usage of gen-colex

This helptext can be found by executing `./gen-colex -h`.

Usage: `./gen-colex n [-h] [-b] [-p] base-graph red-graph blue-graph `

When the -g flag is not set, graphs are read from stdin in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>. When the -g flag is set graphs are read as bitsets from stdin.

The order of the graphs is always the first argument. The blue startgraph argument should always appear after the red startgraph argument. both startgraphs should always appear after the base graph. Otherwise, the order in which the arguments appear does not matter. Be careful not to put an argument immediately after one with an option. E.g. -p#b will not recognise the -b argument.

On.
```
  -b                        blue can select # more edges than red in its turn
  -p                        the starting player in blue if #==2 or red otherwise
  -h                        print help message
  n                         the order of the graph on which the game is played
  base-graph                the base graph with all allowed edges in graph6 format
  red-graph                 the start graph of the red player in graph6 format
  blue-graph                the start graph of the blue player in graph6 format
```

### Examples
`bash Erdos-Game-Generic.sh 7 F~~~w F???? F???? 1 1`
Plays the game using our generator on K_7 with an empty start graph, 1 thread for the solver and red starts (in that order)


`bash Erdos-Game-Cliques.sh 7`
Plays the game using geng on K_7 with an empty start graph, 1 thread for the solver and red starts

`bash Erdos-Game-Generic.sh 7 F~~~o F???? F???? 1 1`
Plays the game using our generator on the colex K_7- with an empty start graph, 1 thread for the solver and red starts

`bash Erdos-Game-Generic.sh 4 C~ C? C? 1 1`
Plays the game using our generator on K_4 with an empty start graph, 1 thread for the solver and red starts

`bash Erdos-Game-Generic.sh 5 D~{ D?? D?? 1 1`
Plays the game using our generator on K_5 with an empty start graph, 1 thread for the solver and red starts

### Verification

The generators can be verified on cliques against eachother with the Compare-Generator.sh file include which compare the output of the two generators with eachother. By running `bash Compare-Generators.sh n` with n the order of the clique on which the game is played the script will return how many non-isomorphic final configurations both implementations generated. Comparing the actual graphs is more difficult as the bitset encoding contains more data (to allow for forbidden edges).

Manually the reader can verify the sizes of the middle layers in both implementations. For comparison the generation of K_4 should return 1-1-2-(3). The last layer size is not printed to gen_error, but to the solver directly. It is also of note that the generator only adds red edges, thus these numbers are not identical to the number of non-isomorphic graphs in the game, where red, blue and uncoloured edges exist. The solver layer sizes should be 3-4-6-4-2-1. These are the number of non-isomorphic graphs after each step of playing the game in reverse.

The sizes for K_5 should be 1-1-2-4-6-(6) on generation and 6-16-37-47-47-30-17-6-2-1 for the solver.

The sizes for K_7 should be 1-1-2-5-10-21-41-65-97-131-148-(148) on generation and 148-1068-4703-14598-31823-56173-74310-84621-74310-58217-35489-19819-8640-3555-1143-366-96-28-7-2-1 for the solver.