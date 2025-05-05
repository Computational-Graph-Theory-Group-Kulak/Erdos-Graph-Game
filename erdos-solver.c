
/**
 * This file implements an algorithm for generating and analyzing edge-colored graphs using the Nauty library.
 * 
 * General Structure:
 * - Includes and Macros (Lines 25-50): Includes necessary libraries and defines macros.
 * - Struct Definitions (Lines 50-85): Defines structures for storing graph information and thread data.
 * - Utility Functions:
 *   - printGraph (Lines 86-119): Prints the graph's edge colors.
 *   - is_clique (Lines 120-134): Checks if a set of vertices forms a clique.
 *   - maxCliques (Lines 135-156): Finds the size of the largest clique.
 *   - determine_fitness_clique (Lines 1357-186): Computes fitness based on clique size.
 *   - determine_fitness_vertex_capture (Lines 187-222): Computes fitness based on vertex capture game.
 *   - determine_fitness_max_degree (Lines 223-250): Computes fitness based on maximum degree.
 *   - readGraph (Lines 251-260): Reads a graph from a graph6 string.
 *   - determine_canonical_labeling (Lines 261-336): Finds the canonical labeling of a graph.
 * - Thread Functions:
 *   - generate_children (Lines 337-450): Generates child graphs from parent graphs in parallel.
 *   - mergesort_children (Lines 451-510): Merges two sorted lists of graphs.
 * - Main Algorithm:
 *   - find_best_game (Lines 511-813): Computes the best game configuration given a set of parent graphs.
 * - Main Function (Lines 814-1027): Parses command-line arguments, reads input graphs, and finds the optimal game configuration.
 */

#include "readGraph/readGraph6.h"
#include "bitset.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include "nauty2_8_6/nauty.h"

 /**
  * A macro that replaces each occurence of indexInEdgeList(vertex_1, vertex_2) by the code after it.
  * This macro computes the colexicographical order of the edges starting from zero.
  * This computes the index of the edge between two vertices as follows: given a vertax a and b with b > a
  * Then the edge index is b*(b-1)/2 (the amount of edges in the sub-clique of order b)+ a
  * for example the edge (1,0) has index 0
  * (2,0) has index 1
  * (2,1) has index 2*1/2+1 = 2
  * (3,2) has index 3*2/2+2 = 5
  * (4,0) has index 4*3/2 = 6
  * 
  */
#define indexInEdgeList(vertex_1, vertex_2) (((vertex_1) < (vertex_2)) ? ((vertex_2) * ((vertex_2) - 1) / 2 + (vertex_1)) : ((vertex_1) * ((vertex_1) - 1) / 2 + (vertex_2)))

#define USAGE "Usage: ./erdos-solver n [-h] [-g] [-t] [-b] [-p] red-graph blue-graph "

/**
 * A struct datatype to store all information of a graph (the variable name graph is already used by Nauty)
 */
struct mygraph
{
    int* fitness;//the fitnesses of this graph
    bitset encoded_canonical; //the bitset storing the encoding of the edge-colored graph
};

/**
 * a struct with all data needed for a graph generation thread
 */
struct thread_data
{
    struct mygraph *parents;//the list of parents
    int parentsnumber;//the number of parents in the list
    bool red_is_playing;//a boolean to store which player is playing
    struct mygraph *children; //the list of generated children
    bitset* children_canonical; //the list of canonical labeling of the children
    int childrennumber; //the number of generated children
    int vertices; //the number of vertices of all graphs
    int fitness_numbers; //the number of fitness values of all graphs
    bitset base_graph; //the bitset with edges that cannot be removed
};


/**
 * a struct with all data needed for a graph merging thread
 * This contains the two sorted lists (canonical and original) and a resulting list
 */
struct mergesort_thread_data
{
    struct mygraph *children_1;
    bitset* children_1_canonical;
    int childrennumber_1;
    struct mygraph *children_2;
    bitset* children_2_canonical;
    int childrennumber_2;
    struct mygraph *children;
    bitset* children_canonical;
    int *childrennumber;
    bool red_is_playing;
    int fitness_numbers;
};


/**
 * a helper function that prints a graph to a specified output in a nice format
 */
void printGraph(FILE *__restrict__ stream, struct mygraph *g, int vertices)
{
    fprintf(stream, "Color 1:\n");
    for (int i = 0; i < vertices; i++)
    {
        for (int j = 0; j < vertices; j++)
        {
            if (i==j){
                fprintf(stream, "n");
            }else{
            
            fprintf(stream, "%d", (g->encoded_canonical & (singleton(2*indexInEdgeList(i, j))))&&(g->encoded_canonical & (singleton (1+2*indexInEdgeList(i, j)))) ? 1 : 0);
            
            }
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "Color 2:\n");
    for (int i = 0; i < vertices; i++)
    {
        for (int j = 0; j < vertices; j++)
        {
            if (i==j){
                fprintf(stream, "n");
            }else{
            
            fprintf(stream, "%d", (g->encoded_canonical & (singleton (2*indexInEdgeList(i, j))))&&((g->encoded_canonical & (singleton (1+2*indexInEdgeList(i, j))))==0)  ? 1 : 0);
            
            }
        }
        fprintf(stream, "\n");
    }
}

/**
 * A function that checks if a given set of vertices (in store) of size b is a clique in a given graph.
 */
bool is_clique(int b, bitset edgeList, int* store)
{

    // Run a loop for all set of edges
    for (int i = 0; i < b; i++) {
        for (int j = i+1; j < b; j++)

            // If any edge is missing
            if ((edgeList & singleton(indexInEdgeList(store[i],store[j]))) == 0){
                //the given set of vertices are not a clique
                return false;
            }
    }
    return true;
}

// Function to find all the sizes
// of maximal cliques
int maxCliques(int i, int l, int n, int* store, bitset edgeList)
{
    // Maximal clique size
    int max_ = 0;

    // Check if any vertices from i+1
    // can be inserted
    for (int j = i+1 ; j < n; j++) {

        // Add the vertex to store
        store[l] = j;

        // If the graph is not a clique of size k then
        // it cannot be a clique by adding another edge
        if (is_clique(l + 1, edgeList, store)) {

            // Update max
            if (max_ < l+1){
                max_ = l+1;
            }
            // Check if another edge can be added
            int temp = maxCliques(j, l + 1, n, store, edgeList);
            if (max_ < temp){
                max_ = temp;
            }
        }
    }
    return max_;
}

int determine_fitness_clique(bitset edgelist, int vertices)
{

    //we split the coloured graph in two non-coloured graphs to compute the largest clique numbers of each colour
    bitset list_red = EMPTY;
    bitset list_blue = EMPTY;
    int edge = 0;
    for(int i = 0;i < vertices;i++){
        for(int j = 0; j < i; j++){
            if(singleton(2*edge) & edgelist){
                if(singleton(2*edge+1) & edgelist){
                    list_red |= singleton(edge);
                }else{
                    list_blue |= singleton(edge);
                }
            }

            edge++;
        }
    }
    //store is a temp memory space needed to check the cliques. it is created once and then reused
    int* store = malloc(sizeof(int)*vertices);
    //using the helper function we compute the largest clique number of red
    int fitness_red = maxCliques(-1, 0, vertices, store, list_red);
    free(store);
    store = malloc(sizeof(int)*vertices);
    //using the helper function we compute the largest clique number of blue
    int fitness_blue = maxCliques(-1, 0, vertices, store, list_blue);
    free(store);
    return 11 * fitness_red - 10 * fitness_blue;
}
/**
 * A function to determine the fitness value of a given game when playing the vertex capture game.
 */
int determine_fitness_vertex_capture(bitset edgelist, int vertices)
{
    
    bitset list_red = EMPTY; // a variable to hold the vertices on which red has won.
    bitset list_blue = EMPTY; // a variable to hold the vertices on which blue has won.
    
    /**
     * A function that for every vertex checks how many edges red coloured and how many blue coloured.
     */
    for(int i = 0;i < vertices;i++){
        int number_A = 0; // the amount of edges of red
        int number_B = 0;//  the amount of edges of blue
        for(int j = 0; j < vertices; j++){
            if(i !=j){ // there is no edge from a vertex to itself
                int edge = indexInEdgeList(i,j);

                if(contains(edgelist,2*edge)){ // the edge is present
                    if(contains(edgelist,2*edge+1)){ // the edge is coloured in colour one
                        number_A +=1;
                    }else{
                        number_B +=1;
                    }
                }
            }

        }
        if(number_A > number_B){ //the current vertex has more incident edges of colour one than of colour two
            add(list_red,i);
        }
        if(number_B > number_A){
            add(list_blue,i);
        }
    }
    return 11*size(list_red)-10*size(list_blue);
}


/**
 * A function to determine the fitness value of a given game when playing the max degree game.
 */
int determine_fitness_max_degree(bitset edgelist, int vertices)
{
    int max_A = 0;// a variable to hold the current maximum degree in the subgraph of red
    int max_B = 0;// a variable to hold the current maximum degree in the subgraph of blue
    
    for(int i = 0;i < vertices;i++){ // check for each vertex
        int number_A = 0; // current number of incident edges in the colour of red
        int number_B = 0;// current number of incident edges in the colour of blue
        for(int j = 0; j < vertices; j++){ //every possible edge
            if(i!=j){ // no edge from a vertex to itself
            int edge = indexInEdgeList(i,j);

            if(contains(edgelist,2*edge)){// if edge is present
                if(contains(edgelist,2*edge+1)){ // if edge is colour of red
                    number_A +=1;
                }else{
                    number_B +=1;
                }
            }
            }

        }
        if(number_A > max_A){
            max_A = number_A;
        }
        if(number_B > max_B){
            max_B = number_B;
        }
    }
    return 11*max_A-10*max_B;
}

/**
 * A helper function to read a graph6 string into a grandparent datatype
 */
int readGraph(const char *graphString, bitset *edgelist)
{
    int vertices = getNumberOfVertices(graphString);

    if (loadGraph(graphString, vertices, edgelist) == -1)
    {
        return 1;
    }
    return 0;
}

/**
 * a function to transform an edge coloured graph into a vertex coloured simple graph
 * every vertex and edge is mapped to a vertex and 4 additional vertices are present for the states of the edges,
 * every vertex is labeled 0 to vertices+ edges +4 in lab
 * ptn contains where each colour class stops
 * We put vertices, edges and the four states each in their own colour class (a 0 denotes the end of a colour class)
 * Vertices are connected to all incident edges and edges to the state that they were
 * 
 */
void generate_expanded_graph(graph *g, int n, int m, int *lab, int*ptn, int vertices, bitset expanded_edgelist){
    EMPTYGRAPH(g, m, n);//
    int edge = 0; // variable to reduce the amount of times the edge-index is computed
    for(int i=0;  i < vertices; i++){ // for every vertex in the graph
            lab[i] = i; // the label of the vertex in the expanded is the label of the vertex in the original graph
            ptn[i] = 1; // All vertices are in the same colour class
        for(int j = 0; j < i; j++){
            lab[vertices+edge] = vertices+edge; // the vertex representing an edge has label #vertices+edge label
            ptn[vertices+edge] = 1; // all edges are in the same colour group
            ADDONEEDGE(g, i, vertices+edge, m); //the vertex representing the vertex incident to the edge 
                                                //and the vertex representing the edge are connected
            ADDONEEDGE(g, j, vertices+edge, m); // similar for the other incident vertex

            /**
             * The next section places an edge between the vertex representing the edge and the corresponding state vertex
             * either not present, present but uncoloured, present coloured red or present coloured blue
             */
            if((expanded_edgelist & (singleton (2*edge)))&&(expanded_edgelist & (singleton (1+2*edge)))){
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2)+1, m);
            }else if((expanded_edgelist & (singleton (2*edge)))&&((expanded_edgelist & (singleton (1+2*edge)))==0)){
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2)+2, m);
            }else if ((expanded_edgelist & (singleton (1+2*edge)))==0){
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2), m);
            } 
            else{
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2)+3, m);
            }


            edge++;// advancing the edge index
        }
    }

    //adding labels to the four state vertices
    lab[vertices+(vertices*(vertices-1)/2)] = vertices+(vertices*(vertices-1)/2);
    lab[vertices+(vertices*(vertices-1)/2)+1] = vertices+(vertices*(vertices-1)/2)+1;
    lab[vertices+(vertices*(vertices-1)/2)+2] = vertices+(vertices*(vertices-1)/2)+2;
    lab[vertices+(vertices*(vertices-1)/2)+3] = vertices+(vertices*(vertices-1)/2)+3;

    //closing the colour classes: one for the vertices, one for the edges and four for each state
    ptn[vertices+(vertices*(vertices-1)/2)] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)+1] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)+2] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)+3] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)-1] = 0;
    ptn[vertices-1] = 0;
}


/**
 * a function to transform a given edgelist of a graph in its canonically labelled form
 * 
 * The function generates the expanded graph with the function generate_expanded_graph
 * and then uses nauty to generate the canonical form
 */
bitset determine_canonical_labeling(bitset expanded_edgelist, int vertices)
{
    int* mapping = malloc(sizeof(int)*vertices); //a variable to hold the new label for each original vertex
    int n = vertices+(vertices*(vertices-1)/2)+4; // the order of the expanded graph is #vertices plus #edges for a complete graph plus four
    int m = SETWORDSNEEDED(n); //m is an internal variable to store the expanded graph in m bitsets
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);//double check if m is set correctly
    graph g[n*m]; //alocate the space for the expanded graph
    graph cg[n*m]; //alocate the space for the resulting canonical graph
    int * lab = malloc(n * sizeof(int)); // a variable to hold the labels of all vertices
    int * ptn = malloc(n * sizeof(int));// a variable to indicate which labels belong to the same colour class
    int orbits[n];  // a variable to store which vertex is the representative of the vertex orbits of each vertex
                    //all vertices with the same representative belong to the same orbit
    static DEFAULTOPTIONS_GRAPH(options); // a standard variable for nauty options with everything set to the default
    options.getcanon = TRUE; // enabling nauty to generate a canonical labeling
    options.defaultptn = FALSE; // enabling the option to provide custom color classes for the vertices
    statsblk stats;// nauty setup

    //generate the expanded graph
    generate_expanded_graph(g,n,m,lab,ptn,vertices,expanded_edgelist);
    
    
    // letting nauty compute the canonical labelling
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);
    

    // transforming the labeling from the output format of nauty to a map from old label to new label
    for (int i = 0; i < vertices; i++) {
            mapping[lab[i]] = i;
    }

    //variable containg the bitset of the canonically labeled graph
    bitset encoded_canonical = EMPTY;

    int edge = 0;//for each edge in order (using colexicographical order)
    for (int i = 0; i < vertices; i++){
        for (int j = 0; j < i; j++){
            // using the vertex mapping, compute where the edges are mapped to
            int mapped_edge = indexInEdgeList(mapping[i],mapping[j]);

            // make the state of the mapped ege the state of the original edge
            if (expanded_edgelist & (singleton (2*edge))) {
                encoded_canonical |= (singleton (2*mapped_edge));
            }
            if(expanded_edgelist & (singleton (1+2*edge))){
                encoded_canonical |= (singleton (2*mapped_edge+1));
            }
            
            edge++;
        }
    }
    /**
     * freeing temporary variables
     */
    free(lab);
    free(ptn);
    free(mapping);
    return encoded_canonical;
}


/**
 * a function to compute the set of optimal children from a set of parent graphs with a fitness
 * this function is run on each thread so its signature is fixed: void return-type and one (void) parameter
 * 
 * The function adds every possible edge to every graph in its list then uses determine_canonical_labeling to find isomorphic copies
 * then the best fitness values (highest is red is playing, lowest if blue is playing) of the two copies is kept for each fitness
 * 
 */
void *generate_children(void *restrict parentslist)
{

    struct thread_data *args = (struct thread_data *)parentslist; // the void parameter is cast to its type: thread_data
    struct mygraph *parents = args->parents; //the parents are extracted from the thread_data
    struct mygraph *children = args->children; // the memory allocated for storing the children
    bitset* children_canonical = args->children_canonical; //the memory allocated for storing the canoncially labelled children
    bitset base_graph = args->base_graph; //the bitset representing the edges that were present at the start of the game (and should not be removed)
    int parentsnumber = args->parentsnumber; //the number of parents in the list

    int vertices = args->vertices;//the number of vertices in each graph

    int fitness_numbers = args->fitness_numbers; //the number of fitness values

    bool red_is_playing = args->red_is_playing;//a variable to store which player is playing


    for (int i = 0; i < parentsnumber; i++) //for each parent
    {
        for (int edge = 0; edge < vertices*(vertices-1)/2; edge++) //for each possible edge
        {
            if (!contains(base_graph,edge) && //the edge cannot be in the base-graph and
            (contains(parents[i].encoded_canonical,2*edge)) && //the edge must be coloured and
            ((red_is_playing && contains(parents[i].encoded_canonical,1+2*edge)) || // Red is playing and the chosen edge is red or
            (!red_is_playing && (contains(parents[i].encoded_canonical,2*edge+1) == 0)))) //Blue is playing and the edge is blue
            {
                //a child is the graph with the selected edge set to uncoloured (removing the coloured bit) but colourable (adding the colourable bit)
                bitset child = parents[i].encoded_canonical;
                removeElement(child,2*edge);
                add(child,2*edge+1);
                //finding the canonical labeling
                bitset encoded_canonical_new = determine_canonical_labeling(child, vertices);

                /**
                 * The next section is a binary search
                 */
                    int search_index = args->childrennumber/2;
                    int begin = 0;
                    int end = args->childrennumber-1;
                    while ((begin <= end) && (children_canonical[search_index] != encoded_canonical_new))
                    {
                        if (children_canonical[search_index] < encoded_canonical_new)
                        {
                            begin = search_index + 1;
                        }
                        else
                        {
                            end = search_index - 1;
                        }
                        search_index = (begin + end) / 2;
                    }
                    //The search has stopped and we have found a graph with identical canonical form
                    if  (search_index < args->childrennumber && children_canonical[search_index] == encoded_canonical_new)
                        {
                            //compare every fitness function
                            for(int fitness_i =0; fitness_i < fitness_numbers; fitness_i++){
                            //if red was playing we select the highest fitness value
                            if (red_is_playing)
                            {
                                if (parents[i].fitness[fitness_i] > children[search_index].fitness[fitness_i])
                                {
                                    children[search_index].fitness[fitness_i] = parents[i].fitness[fitness_i];
                                    
                                }
                            }
                            //if blue was playing select the lowest fitness value
                            else
                            {
                                if (parents[i].fitness[fitness_i] < children[search_index].fitness[fitness_i])
                                {
                                    children[search_index].fitness[fitness_i] = parents[i].fitness[fitness_i];
                                }
                            }
                            }
                        }
                    //If the search was completed and no isomorphic graph was found at the place where it should be
                    //add in the new graph.
                    else
                    {
                        // if the search ended in the middle of the list we have to position the new graph correctly.
                        if(args->childrennumber != 0 && children_canonical[search_index] < encoded_canonical_new){
                            search_index +=1;
                        }

                        // moving all graphs after the current one to the right to make space
                        for(int j = args->childrennumber; j > search_index; j--){
                            children[j] = children[j-1];
                            children_canonical[j] = children_canonical[j-1];
                        }
                        //adding the new graph in both the canonical list as well as the list with its original form
                        children_canonical[search_index] = encoded_canonical_new;
                        children[search_index].fitness = malloc(sizeof(int)*fitness_numbers);
                        for(int fitness_i =0; fitness_i < fitness_numbers;fitness_i++){
                        children[search_index].fitness[fitness_i] = parents[i].fitness[fitness_i];
                        }
                        //by storing the original labelling we keep the orientation identical to the base-graph
                        children[search_index].encoded_canonical =child;
                        args->childrennumber++;
                    }
            }
        }
    }

    return NULL;
}


/**
 * A function to take two lists of graphs with their canonical labelling lexicographically ordered and merge them into
 * a single list with the correct fitness value
 * this function is run on each thread so its signature is fixed: void return-type and one (void) parameter
 */
void *mergesort_children(void *childrenlist){
    struct mergesort_thread_data *args = (struct mergesort_thread_data *)childrenlist;//casting the argument to the correct type
    struct mygraph *children_1 = args->children_1; //the first sorted list of graphs in their original labelling
    bitset* children_1_canonical = args->children_1_canonical; //the first sorted list of graph represnetations in their canonical labelling
    struct mygraph *children_2 = args->children_2; //the second sorted list of graphs in their original labelling
    bitset* children_2_canonical = args->children_2_canonical; //the second sorted list of graph represnetations in their canonical labelling
    struct mygraph *children = args->children; //the resulting sorted list of graphs in their original labelling
    bitset* children_canonical = args->children_canonical;//the resulting sorted list of graph represnetations in their canonical labelling
    bool red_is_playing = args->red_is_playing;
    int* childrennumber = args->childrennumber; //the number of resulting graphs

    int fitness_numbers = args->fitness_numbers;//the number of fitness values
    
    int childrennumber_1 = args->childrennumber_1; //the number of graphs in the first list
    int childrennumber_2 = args->childrennumber_2; //the number of graphs in the second list

    /**
     * The next section is a merge step as in a mergesort algorithm
     */
    int pointer_1 = 0; //the number of graphs added from the first list
    int pointer_2 = 0;//the number of graphs added from the second list
    while(pointer_1 < childrennumber_1 && pointer_2 < childrennumber_2){
        if(children_1_canonical[pointer_1] < children_2_canonical[pointer_2]){
            children[*childrennumber] = children_1[pointer_1];
            children_canonical[*childrennumber]= children_1_canonical[pointer_1];
            pointer_1++;
        }else if(children_1_canonical[pointer_1] > children_2_canonical[pointer_2]){
            children[*childrennumber] = children_2[pointer_2];
            children_canonical[*childrennumber]= children_2_canonical[pointer_2];
            pointer_2++;
        }else{
            children[*childrennumber] = children_1[pointer_1];
            children_canonical[*childrennumber]= children_1_canonical[pointer_1];
            for(int fitness_i =0; fitness_i <fitness_numbers;fitness_i++){
            if(red_is_playing){
                if(children_1[pointer_1].fitness[fitness_i] > children_2[pointer_2].fitness[fitness_i]){
                    children[*childrennumber].fitness[fitness_i] = children_1[pointer_1].fitness[fitness_i];
                }
                else{
                    children[*childrennumber].fitness[fitness_i] = children_2[pointer_2].fitness[fitness_i];
                    
                }
            }else{
                if(children_1[pointer_1].fitness[fitness_i] > children_2[pointer_2].fitness[fitness_i]){
                    children[*childrennumber].fitness[fitness_i] = children_2[pointer_2].fitness[fitness_i];
                }else{
                    children[*childrennumber].fitness[fitness_i] = children_1[pointer_1].fitness[fitness_i];
                }
            }
            }
            free(children_2[pointer_2].fitness);
            pointer_1++;
            pointer_2++;
        }
        (*childrennumber)++;
    }
    while(pointer_1 < childrennumber_1){
        children[*childrennumber] = children_1[pointer_1];
        children_canonical[*childrennumber]= children_1_canonical[pointer_1];
        pointer_1++;
        (*childrennumber)++;
    }
    while(pointer_2 < childrennumber_2){
        children[*childrennumber] = children_2[pointer_2];
        children_canonical[*childrennumber]= children_2_canonical[pointer_2];
        pointer_2++;
        (*childrennumber)++;
    }
    return NULL;
}
/**
 * The main function to compute what the score of the best game would be given a set of parents
 * @param parents : a pointer to the list of parents the game starts with with their fitness set to the score an optimal game starting in the graph would lead to
 * @param parentsnumber : the number of parents in the list of parents
 * @param red_is_playing : a boolean expressing if the first player is currently playing
 * @param max_threads : the maximum number of threads that can be started to speed-up computation
 * @param children_in_single_thread : the amount of children that should minimally be generated on a single thread
 * @param vertices : the number of vertices in the game
 * @param base_graph : an edgelist containing all edges that were present in the starting configuration and cannot be removed
 * @param fitness_numbers : the number of fitness values per graph
 * 
 * 
 * @return : the starting graph of the game with the fitness value set to the score an optimal game would lead to
 * 
 * This function checks each time if there are enough parents to split over multiple threads, 
 * in that case it uses generate_children and mergesort_children to construct the list of children over multiple threads
 * otherwise is performs the same logic as generate_children on this thread (no merging required)
 *
 */
struct mygraph find_best_game(struct mygraph *parents, int parentsnumber, bool red_is_playing, unsigned long max_threads,int children_in_single_thread, int vertices, bitset base_graph, int fitness_numbers)
{
    //first step is computing how many threads will be used. This is to make a trade-off between the overhead cost and performance gain
    // We use the number of children generated as an approximation for the amount of work to be done. However it must be noted that isomorphic copies take less time than new graphs
    int children_per_graph = (size(parents[0].encoded_canonical) / 3 + 1);
    int threads = parentsnumber * children_per_graph / children_in_single_thread > max_threads ? max_threads : parentsnumber * children_per_graph / children_in_single_thread;
    int childrennumber = 0; //setting the number of children generated to 0
    struct mygraph *children = malloc(sizeof(struct mygraph) * parentsnumber * children_per_graph);; //alocating a variable for the generated children
    bitset* children_canonical = malloc(sizeof(bitset) * parentsnumber * children_per_graph);
    if (threads > 1) //if it was usefull to multithread based on the given max number of threads and the minimum number of graphs per thread
    {
        int parents_per_thread = parentsnumber / threads; //determine the amount of parents per thread given the selected number of threads.
        fprintf(stderr, "parents_per_thread: %d\n", parents_per_thread);
        struct thread_data *args = malloc(sizeof(struct thread_data) * threads);//allocate space for the necessary thread-data
        pthread_t *thread = malloc(sizeof(pthread_t) * threads);//allocate space for the pointers to the pthreads
        for (int t = 0; t < threads; t++)//for each thread
        {
            if (t == threads - 1) //the last thread has to take the remainder of the parents
            {
                //allocate space for the necessary parent graphs
                args[t].parents = parents + t * parents_per_thread;
                //set the correct number of parents
                args[t].parentsnumber = parents_per_thread + parentsnumber % threads;
                //allocate a worst-case expectation of the space needed to store the non-isomorphic children
                args[t].children = children + t * parents_per_thread* children_per_graph;
                args[t].children_canonical = children_canonical + t * parents_per_thread* children_per_graph;
            }
            else //every other thread has an equal number of parents
            {
                //allocate space for the necessary parent graphs
                args[t].parents = parents + t * parents_per_thread;
                //set the correct number of parents
                args[t].parentsnumber = parents_per_thread;
                //allocate a worst-case expectation of the space needed to store the non-isomorphic children
                args[t].children = children + t * parents_per_thread*children_per_graph;
                args[t].children_canonical = children_canonical + t * parents_per_thread*children_per_graph;
            }
            //store basic information about the game for each thread
            args[t].red_is_playing = red_is_playing;
            args[t].childrennumber = 0;
            args[t].vertices = vertices;
            args[t].fitness_numbers = fitness_numbers;
            args[t].base_graph = base_graph;
            //start the threads with the generate_children function and the data prepared above
            pthread_create(&thread[t], NULL, generate_children, &args[t]);
        }
        //wait until every thread has finished
        for (int t = 0; t < threads; t++)
        {
            pthread_join(thread[t], NULL);
        }
        //the space to store the threads themselves can be freed (not their data)
        free(thread);
        //the space to store the previous layer can be freed
        for(int p = 0;p<parentsnumber;p++){
            free(parents[p].fitness);
        }
        free(parents);
        //We use one thread for every pair of threads to merge the t independent sets of sorted children
        int new_threads = (threads+1) / 2;

        bool even = ((threads % 2) == 0);
        //allocate space for the threads to merge all the independent lists
        //oldargs and newargs will be used to refer to the data to merge and the data after merging, as the merging will happen iteratively
        struct mergesort_thread_data *oldargs = malloc(sizeof(struct mergesort_thread_data) * new_threads);
        thread = malloc(sizeof(pthread_t) * new_threads);
        /**
         * prepare all merge data by transforming it from generation-thread data to merge-thread data
         */
        for(int t = 0; t < new_threads; t++){
            oldargs[t].red_is_playing = red_is_playing;
            //all threads except for the last one in the case that we have an odd number of threads merge two of the previous lists
            if(even || t != new_threads-1){
            //copying the graphs from the previous datastructure to the new.
            oldargs[t].children_1 = malloc(sizeof(struct mygraph) * args[2*t].childrennumber);
            oldargs[t].children_2 = malloc(sizeof(struct mygraph) * args[2*t+1].childrennumber);
            oldargs[t].childrennumber_1 = args[2*t].childrennumber;
            oldargs[t].childrennumber_2 = args[2*t+1].childrennumber;
            oldargs[t].children_1_canonical = malloc(sizeof(bitset) * args[2*t].childrennumber);
            oldargs[t].children_2_canonical = malloc(sizeof(bitset) * args[2*t+1].childrennumber);
            for(int i = 0; i < args[2*t].childrennumber; i++){
                oldargs[t].children_1[i] = args[2*t].children[i];
                oldargs[t].children_1_canonical[i] = args[2*t].children_canonical[i];
            }
            for(int i = 0; i < args[2*t+1].childrennumber; i++){
                oldargs[t].children_2[i] = args[2*t+1].children[i];
                oldargs[t].children_2_canonical[i] = args[2*t+1].children_canonical[i];
            }

            //allocating a worst-case space for the result of the merge
            oldargs[t].children = malloc(sizeof(struct mygraph) * (args[2*t].childrennumber + args[2*t+1].childrennumber));
            oldargs[t].children_canonical = malloc(sizeof(bitset) * (((args[2*t].childrennumber + args[2*t+1].childrennumber) < (parentsnumber * children_per_graph)) ? (args[2*t].childrennumber + args[2*t+1].childrennumber) : (parentsnumber * children_per_graph)));
            oldargs[t].childrennumber = malloc(sizeof(int));
            *(oldargs[t].childrennumber) = 0;
            oldargs[t].fitness_numbers = fitness_numbers;
            
            //running the merge thread
            pthread_create(&thread[t], NULL, mergesort_children, &oldargs[t]);
            }else //if the number of threads was odd, then the last thread should not be merged, but the data must be put in the correct format
            {
                //setting the resulting list to the original list as there is no list to merge with.
                oldargs[t].children = malloc(sizeof(struct mygraph) * (args[2*t].childrennumber));
                oldargs[t].children_canonical = malloc(sizeof(bitset) * (args[2*t].childrennumber));
                oldargs[t].childrennumber = malloc(sizeof(int));
                *(oldargs[t].childrennumber) = args[2*t].childrennumber;
                oldargs[t].childrennumber_1 = args[2*t].childrennumber;
                oldargs[t].childrennumber_2 = args[2*t].childrennumber;
                for(int i = 0; i < args[2*t].childrennumber; i++){
                    oldargs[t].children[i] = args[2*t].children[i];
                    oldargs[t].children_canonical[i] = args[2*t].children_canonical[i];
                }
            }
        }
        /**
         * waiting for all threads (except for the last when the number was odd) to finish
         */
        for (int t = 0; t < new_threads; t++)
           {
            if(even || t != new_threads-1){
               pthread_join(thread[t], NULL);
            }
           }
        free(thread);
        for(int i = 0;i<childrennumber;i++){
            free(children[i].fitness);
        }
        free(children);

        //freeing the allocated space for the input
        for (int t = 0; t < new_threads; t++){
            if(even || t != new_threads-1){
                free(oldargs[t].children_1);
                free(oldargs[t].children_2);
                free(oldargs[t].children_1_canonical);
                free(oldargs[t].children_2_canonical);
                }
           }
        

           /**
            * repeat the process above, until only a single list is constructed
            */
        while (new_threads > 1)
        {
            even = ((new_threads % 2) == 0);
            new_threads = (new_threads+1) / 2;
            struct mergesort_thread_data *newargs = malloc(sizeof(struct mergesort_thread_data) * new_threads);
            thread = malloc(sizeof(pthread_t) * new_threads);
            for (int t = 0; t < new_threads; t++)
            {
                newargs[t].red_is_playing = red_is_playing;
                if(even || t != new_threads-1){
                    newargs[t].children_1 = malloc(*oldargs[2 * t].childrennumber * sizeof(struct mygraph));
                    newargs[t].children_2 = malloc(*oldargs[2 * t + 1].childrennumber * sizeof(struct mygraph));
                    newargs[t].children_1_canonical = malloc(*oldargs[2 * t].childrennumber * sizeof(bitset));
                    newargs[t].children_2_canonical = malloc(*oldargs[2 * t + 1].childrennumber * sizeof(bitset));
                    for(int i = 0; i < *oldargs[2 * t].childrennumber; i++){
                        newargs[t].children_1[i] = oldargs[2 * t].children[i];
                        newargs[t].children_1_canonical[i] = oldargs[2 * t].children_canonical[i];
                    }
                    free(oldargs[2 * t].children);
                    free(oldargs[2 * t].children_canonical);
                    for(int i = 0; i < *oldargs[2 * t + 1].childrennumber; i++){
                        newargs[t].children_2[i] = oldargs[2 * t + 1].children[i];
                        newargs[t].children_2_canonical[i] = oldargs[2 * t + 1].children_canonical[i];
                    }
                    free(oldargs[2 * t + 1].children);
                    free(oldargs[2 * t + 1].children_canonical);
               
                    newargs[t].children = malloc(sizeof(struct mygraph) * (*oldargs[2 * t].childrennumber + *oldargs[2 * t + 1].childrennumber));
                    newargs[t].children_canonical = malloc(sizeof(bitset) * (*oldargs[2 * t].childrennumber + *oldargs[2 * t + 1].childrennumber));
                    newargs[t].childrennumber = malloc(sizeof(int));
                    *(newargs[t].childrennumber) = 0;
                    newargs[t].childrennumber_1 = *oldargs[2 * t].childrennumber;
                    newargs[t].childrennumber_2 = *oldargs[2 * t + 1].childrennumber;
                    free(oldargs[2 * t].childrennumber);
                    free(oldargs[2 * t + 1].childrennumber);
                    newargs[t].fitness_numbers = fitness_numbers;
                    pthread_create(&thread[t], NULL, mergesort_children, &newargs[t]);
                }else{
                    
                    newargs[t].children = malloc(sizeof(struct mygraph) * (*oldargs[2*t].childrennumber));
                    newargs[t].children_canonical = malloc(sizeof(bitset) * (*oldargs[2*t].childrennumber));
                    newargs[t].childrennumber = malloc(sizeof(int));
                    *(newargs[t].childrennumber) = *oldargs[2*t].childrennumber;
                
                    for(int i = 0; i < *oldargs[2*t].childrennumber; i++){
                        newargs[t].children[i] = oldargs[2*t].children[i];
                        newargs[t].children_canonical[i] = oldargs[2*t].children_canonical[i];
                    }
                    free(oldargs[2 * t].children);
                    free(oldargs[2 * t].children_canonical);
                    free(oldargs[2 * t].childrennumber);
                    }
                }
            free(oldargs);
            for (int t = 0; t < new_threads; t++)
            {
                if(even || t != new_threads-1){
                    pthread_join(thread[t], NULL);
                }
            }
            for (int t = 0; t < new_threads; t++)
            {
                if(even || t != new_threads-1){
                    free(newargs[t].children_1);
                    free(newargs[t].children_2);
                    free(newargs[t].children_1_canonical);
                    free(newargs[t].children_2_canonical);
                }
            }
            free(thread);
            oldargs = newargs;
        
        }

        //when a single list is constructed we extract the sorted list from the arguments
        children = oldargs[0].children;
        childrennumber = *oldargs[0].childrennumber;
       
        //we can free all other datatypes
        free(oldargs[0].children_canonical);
        free(oldargs[0].childrennumber);
        free(oldargs);
        free(args);
        
        
    }
    else //If this step is run in a single thread
    {
        for (int i = 0; i < parentsnumber; i++)//for each existing parent graph
        {
            for (int edge = 0; edge < vertices*(vertices-1)/2; edge++)//for each edge
            {
                
                if (!contains(base_graph,edge) && //the edge cannot be in the base-graph and
                (contains(parents[i].encoded_canonical,2*edge)) && //the edge must be coloured and
                ((red_is_playing && contains(parents[i].encoded_canonical,1+2*edge)) || // Red is playing and the chosen edge is red or
                (!red_is_playing && (contains(parents[i].encoded_canonical,2*edge+1) == 0)))) //Blue is playing and the edge is blue
                    {
                    //a child is the graph with the selected edge set to uncoloured (removing the coloured bit) but colourable (adding the colourable bit)
                    bitset child = parents[i].encoded_canonical;
                    removeElement(child,2*edge);
                    add(child,2*edge+1);

                    //finding the canonical labeling
                    bitset encoded_canonical = determine_canonical_labeling(child, vertices);

                    /**
                     * The next section is a binary search
                     */
                    int search_index = childrennumber/2;
                    int begin = 0;
                    int end = childrennumber-1;
                    while ((begin <= end) && (children_canonical[search_index] != encoded_canonical))
                    {
                        if (children_canonical[search_index] < encoded_canonical)
                        {
                            begin = search_index + 1;
                        }
                        else
                        {
                            end = search_index - 1;
                        }
                        search_index = (begin + end) / 2;
                    }

                    //The search has stopped and we have found a graph with identical canonical form
                    if  (search_index < childrennumber && children_canonical[search_index] == encoded_canonical)
                        {
                            //compare every fitness function
                            for(int fitness_i =0; fitness_i <fitness_numbers;fitness_i++){
                            //if red was playing we select the highest fitness value
                            if (red_is_playing)
                            {
                                if (parents[i].fitness[fitness_i] > children[search_index].fitness[fitness_i])
                                {
                                    children[search_index].fitness[fitness_i] = parents[i].fitness[fitness_i];
                                    
                                }
                            }
                            //if blue was playing select the lowest fitness value
                            else
                            {
                                if (parents[i].fitness[fitness_i] < children[search_index].fitness[fitness_i])
                                {
                                    children[search_index].fitness[fitness_i] = parents[i].fitness[fitness_i];
                                }
                            }
                            }
                        }
                    //If the search was completed and no isomorphic graph was found at the place where it should be
                    //add in the new graph.
                    else
                    {
                        // if the search ended in the middle of the list we have to position the new graph correctly.
                        if(childrennumber != 0 && children_canonical[search_index] < encoded_canonical){
                            search_index +=1;
                        }

                        // moving all graphs after the current one to the right to make space
                        for(int j = childrennumber; j > search_index; j--){
                            children[j] = children[j-1];
                            children_canonical[j] = children_canonical[j-1];
                        }
                        //adding the new graph in both the canonical list as well as the list with its original form
                        children_canonical[search_index] = encoded_canonical;
                        children[search_index].fitness = malloc(sizeof(int)*fitness_numbers);
                        for(int fitness_i =0; fitness_i <fitness_numbers;fitness_i++){
                        children[search_index].fitness[fitness_i] = parents[i].fitness[fitness_i];
                        }
                        //by storing the original labelling we keep the orientation identical to the base-graph
                        children[search_index].encoded_canonical =child;
                        childrennumber++;
                    }
        
                }
            }
            //free the fitness array of the parent
            free(parents[i].fitness);
        }
        //free the memory for the list of previous graphs
        free(parents);
    }
    //free the list of canonical labellings
    free(children_canonical);

    //print layer sizes and winning percentage for debugging purposes
    fprintf(stderr, "childrennumber: %d\n", childrennumber);
    for (int fitness_i =0; fitness_i <fitness_numbers; fitness_i ++){
        int percentage_winning = 0;
        for(int winning_graph =0; winning_graph < childrennumber;winning_graph++){
            if(children[winning_graph].fitness[fitness_i] > 10){
                percentage_winning ++;
            }
        }
        fprintf(stderr, "percentage winning %i: %d out of %d\n",fitness_i, percentage_winning, childrennumber);
    }
    //if we have converged to a single non-isomorphic child then no more moves matter
    // and we can print the fitness of this option as the score of the optimal game
    if (childrennumber == 1)
    {
        //extract the single graph from the list and free the list
        struct mygraph return_value = children[0];
        free(children);
        return return_value;
    }
    else //if we still have more options we play the game again
    {   
        return find_best_game(children, childrennumber, !red_is_playing, max_threads,children_in_single_thread,vertices,base_graph,fitness_numbers);
    }
}

/**
 * The main function taking multiple command line arguments
 */
int main(int argc, char **argv)
{

    unsigned long vertices = strtol(argv[1], NULL, 10); //the first argument is the number of vertices in the games
    unsigned long threads = 1; //the second argument is the maximum amount of threads that can be used
    bitset base_graph_grandparent = EMPTY; //The third argument is a graph6 string with the red subgraph that was present before the players began
    bitset second_base_graph_grandparent = EMPTY; //The fourth argument is a graph6 string with the blue subgraph that was present before the players began
    //The fifth argument switches between the encoded format for the given end-configurations or the graphs6 string of the red subgraphs
    bool uses_generator = false;
    //the sixth argument is if there was a bias between how many edges blue colored in its turn and how many red colored
    unsigned long bias = 0; //the sixth argument is if there was a bias between how many edges blue colored in its turn and how many red colored

    bool first_player = true; //an optional sixth argument is which player started the game, if no first player is given red is taken as default
    bool selected_first_player = false;
    
    bool red_base_graph_given = false;
    
    for(int i =2; i < argc;i++){
        if(argv[i][0] == '-'){

            if(argv[i][1] == 'h'){
                fprintf(stderr,USAGE);
                exit(-1);
            }

            if(argv[i][1] == 'g'){
                uses_generator = true;
            }

            if(argv[i][1] == 't'){
                threads = strtol(&argv[i][2], NULL, 10);
            }

            if(argv[i][1] == 'b'){
                bias += strtol(&argv[i][2], NULL, 10);
            }

            if(argv[i][1] == 'p'){
                selected_first_player = true;
                if (strtol(&argv[i][2], NULL, 10)==2){
                    first_player = false;
                }
            }
        }else{
        if(!red_base_graph_given){
            //To read the starting graphs in graph6 format they must end with an endofline '\n' character. 
            //This is not parsed correctly, thus these are added in
            char *base_graph_str_complete = malloc(strlen(argv[i]) + 2);
            strcpy(base_graph_str_complete, argv[i]);
            strcat(base_graph_str_complete, "\n");
            
            //then the graph6 string is read
            readGraph(base_graph_str_complete,&base_graph_grandparent);
            //and the memory for the extended notation freed
            free(base_graph_str_complete);
            red_base_graph_given = true;

        }else{
            char *second_base_graph_str_complete = malloc(strlen(argv[i]) + 2);
            strcpy(second_base_graph_str_complete, argv[i]);
            strcat(second_base_graph_str_complete, "\n");

            readGraph(second_base_graph_str_complete,&second_base_graph_grandparent);

            free(second_base_graph_str_complete);
        }
        }

    }

    //both edgelists are printed for debugging purposes
    fprintf(stderr, "base_graph.edgelist: ");
    for(int i = 0; i < vertices*(vertices-1); i++){
            fprintf(stderr, "%d", (base_graph_grandparent & (singleton (i))) ? 1 : 0);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "second_base_graph.edgelist: ");
    for(int i = 0; i < vertices*(vertices-1); i++){
            fprintf(stderr, "%d", (second_base_graph_grandparent & (singleton (i))) ? 1 : 0);
    }
    fprintf(stderr, "\n");

    
    
    
    //whether this value is supplied is printed for debugging purposes
    fprintf(stderr,"selected first player: %s\n",(selected_first_player?"true":"false"));

    
    if (!selected_first_player){

        //if no player is supplied and there is a non-empty red starting graph, blue is taken as first player
        if(size(base_graph_grandparent)>0){
            //if however there is also a non-empty blue starting graph an error is thrown and no assumption is made on the starting player
            if(size(second_base_graph_grandparent) != 0){
                fprintf(stderr,"missing first player\n");
                exit(-1);
            }
            first_player = false;
                
        }
    }
    
    
    //timing is started here
    clock_t start = clock();

    //  Start looping over lines of stdin.
    //the first number supplied should be an overestimation of the number of graphs that will follow
    char *graphString = NULL;
    size_t size;
    getline(&graphString, &size, stdin);
    long graphnumber_overestimation = strtol(graphString,NULL,10);

    //space is allocated according to this provided number.
    bitset *grandparents_overestimation = malloc(sizeof(bitset) * graphnumber_overestimation);
    int number_of_graphs= 0;
    fprintf(stderr, "graphnumber: %ld\n", graphnumber_overestimation);
    //read all graphs supplied
    while (getline(&graphString, &size, stdin) != -1)
    {

        if(!uses_generator){    //if the -b flag was not set graphs are expected to be in graph6 format
                                //these graphs contain the red subgraph and it is assumed the other edges in the complete graphs were blue
            int status = readGraph(graphString, &grandparents_overestimation[number_of_graphs]);
            if(status == 0){ //if no error occured
                number_of_graphs++;
            }
        }
        else{ //if the -b flag was set the graphs are assumed to be bitset (read as unsigned longs) with two bits per edge as explained in the encoding
            grandparents_overestimation[number_of_graphs] = strtol(graphString, NULL, 10);
            number_of_graphs++;
        }
    }
    //space is allocated to correctly store all graphs after they were counted
    struct mygraph *parents = malloc(sizeof(struct mygraph) * number_of_graphs);
    int distinct_graphs = 0;

    //If the graph were supplied in graph6 format, they are expanded to their bitset encoding and all non-red edges are coloured blue.
    //if they were already in bitset format
    int edges =0;
    if(!uses_generator){
        edges = vertices*(vertices-1)/2; //the number of edges in the complete graph
    }
    for (int graph = 0; graph < number_of_graphs; graph++)//for each graph
    {
        if(!uses_generator){
            //start with a graph with all edges present and colored red
            bitset temp_encoded_canonical = (singleton(vertices*(vertices-1)) - 1);
            
            //for each edge
            int edge = 0;
            for(int i= 0; i < vertices ;i++){
                for(int j=0; j< i;j++){
                
                    if((grandparents_overestimation[graph] & singleton(edge)) == 0){
                        //if it was not present in the input graph, color it blue
                        removeElement(temp_encoded_canonical, 2*edge+1);
                    }
                    edge ++;
                
                }
            }
            //canonically label the graphs to detect isomorphic copies
            temp_encoded_canonical = determine_canonical_labeling(temp_encoded_canonical, vertices);
            

            //search for the graph in the list with distinct graphs
            int temp_graph = 0;
            while(temp_graph < distinct_graphs && parents[temp_graph].encoded_canonical < temp_encoded_canonical){
                temp_graph++;
            }
            //if it is not yet in the list, because all graphs have a smaller encoding or there is a gap in the list were this encoding should fit
            if(temp_graph == distinct_graphs || (temp_graph < distinct_graphs && parents[temp_graph].encoded_canonical > temp_encoded_canonical)){
                //move all graphs to the right of the new graph
                for (int moving_graph = distinct_graphs; moving_graph > temp_graph; moving_graph--){
                    parents[moving_graph] = parents[moving_graph-1];
                }
                //insert the new graph and compute all fitness values
                parents[temp_graph].encoded_canonical = temp_encoded_canonical;
                parents[temp_graph].fitness = malloc(sizeof(int)*3);
                parents[temp_graph].fitness[0]= determine_fitness_clique(temp_encoded_canonical,vertices);
                parents[temp_graph].fitness[1]= determine_fitness_vertex_capture(temp_encoded_canonical,vertices);
                parents[temp_graph].fitness[2]= determine_fitness_max_degree(temp_encoded_canonical,vertices);
        
                distinct_graphs++;
            }
        }else{
            //if we generated the graphs, then we assume that graphs were not isomorphic and we maintain their labeling
            //this is to allow reference to the provided red and blue starting graphs
            //first we sort them with the given encoding
            int temp_graph = 0;
            while(temp_graph < distinct_graphs && parents[temp_graph].encoded_canonical < grandparents_overestimation[graph]){
                temp_graph++;
            }
            if(temp_graph == distinct_graphs || (temp_graph < distinct_graphs && parents[temp_graph].encoded_canonical > grandparents_overestimation[graph])){
                for (int moving_graph = distinct_graphs; moving_graph > temp_graph; moving_graph--){
                    parents[moving_graph] = parents[moving_graph-1];
                }
                //fitness values are computed and added
                parents[temp_graph].encoded_canonical = grandparents_overestimation[graph];
                parents[temp_graph].fitness = malloc(sizeof(int)*3);
                parents[temp_graph].fitness[0]= determine_fitness_clique(grandparents_overestimation[graph],vertices);
                parents[temp_graph].fitness[1]= determine_fitness_vertex_capture(grandparents_overestimation[graph],vertices);
                parents[temp_graph].fitness[2]= determine_fitness_max_degree(grandparents_overestimation[graph],vertices);
           
                distinct_graphs++;
            }
        }
    }

    free(grandparents_overestimation); // free the space we allocated at the beginning.

    //if the games were generated using the supplied generator we computed the number of edges that were actually played
    //first we count the edges that are present
    if(uses_generator){
        edges = 0;
        for(int edge = 0; edge < (vertices*(vertices-1))/2;edge++){
            if(parents[0].encoded_canonical & singleton(2*edge)){
                edges++;
            }
        }
    }

    //then we initialize a variable to store how many edges blue would have colored if he was the last player to choose edges
    int leftover_edges =  0;
    //initializing a variable to track which player made the last move
    bool last_player_is_red = first_player;
    if(first_player){
        //the edges that were actually used in the game are the edges that could have been used (counted in the variable edges)
        //minus the edges in the two base-graphs

        //Then the number of edges in the last move if red started are the number of edges actually played minus one modulo 2+bias
        //where bias is the difference in edges played in each turn between red and blue
        //the minus one is becasue the first edge on the turn is colored by red
        if((edges-size(base_graph_grandparent)-size(second_base_graph_grandparent)-1)%(2+bias)!=0){

            last_player_is_red = false;
            
            leftover_edges = (edges-size(base_graph_grandparent)-size(second_base_graph_grandparent)-1)%(2+bias);

        }
    }else{

        //Then the number of edges in the last move if blue started are the number of edges actually played modulo 2+bias
        //where bias is the difference in edges played in each turn between red and blue
        if((edges-size(base_graph_grandparent)-size(second_base_graph_grandparent))%(2+bias)==0){
            
            //change the  value of last_player_is_red from false to true
            last_player_is_red = true;
            
        
        }
        else{


            leftover_edges = (edges-size(base_graph_grandparent)-size(second_base_graph_grandparent))%(2+bias);
        }
    }
    fprintf(stderr,"last player was player: %s\n",(last_player_is_red?"one":"two"));

    //for the solver it is not important which color the edge has, just that it was in the start-configuration and cannot be altered
    //therefore a new variable is created to store all edges present in the base (starting) graph
    bitset combined_base_graph = base_graph_grandparent | second_base_graph_grandparent;
    fprintf(stderr, "combined_base_graph.edgelist: ");
    for(int i = 0; i < vertices*(vertices-1); i++){
            fprintf(stderr, "%d", (combined_base_graph & (singleton (i))) ? 1 : 0);
    }
    fprintf(stderr, "\n");

    //call the function on the initialized game
    struct mygraph winner = find_best_game( parents, distinct_graphs, last_player_is_red, threads, 3000000, vertices, combined_base_graph,3);
    
    //print all information about the winning graph
    fprintf(stdout, "The best graph is:\n");
    printGraph(stdout,&winner,vertices);
    fprintf(stdout, "With fitness clique: %d\n", winner.fitness[0]);
    fprintf(stdout, "With fitness vertex capture: %d\n", winner.fitness[1]);
    fprintf(stdout, "With fitness max degree: %d\n", winner.fitness[2]);

    //stop the timing
    clock_t end = clock();

    //compute the total running time of the process(es)
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    //print the running time
    fprintf(stderr, "\rChecked in %f seconds:", time_spent);

    //free last variables

    free(winner.fitness);

    return 0;
}
