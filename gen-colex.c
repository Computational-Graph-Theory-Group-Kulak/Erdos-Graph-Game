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
#include "writeGraph/writeGraph6.h"


/**
  * A macro that replaces each occurrence of indexInEdgeList(vertex_1, vertex_2) by the code after it.
  * This macro computes the colexicographical order of the edges starting from zero.
  * This computes the index of the edge between two vertices as follows: given a vertex a and b with b > a
  * Then the edge index is b*(b-1)/2 (the amount of edges in the sub-clique of order b)+ a
  * for example the edge (1,0) has index 0
  * (2,0) has index 1
  * (2,1) has index 2*1/2+1 = 2
  * (3,2) has index 3*2/2+2 = 5
  * (4,0) has index 4*3/2 = 6
  * 
  */
#define indexInEdgeList(vertex_1, vertex_2) (((vertex_1) < (vertex_2)) ? ((vertex_2) * ((vertex_2) - 1) / 2 + (vertex_1)) : ((vertex_1) * ((vertex_1) - 1) / 2 + (vertex_2)))

#define USAGE "Usage: ./gen-colex n [-h] [-b] [-p] base-graph red-graph blue-graph "

/**
 * A struct datatype to store all information of a graph (the variable name graph is already used by Nauty)
 */
struct mygraph
{
    int numberOfVertices;
    bitset edgelist;
    bitset possible_edges;
    int last_edge;

};

/**
 * a helper function to read a string and convert it into a mygraph struct
 */
int readGraph(const char *graphString, bitset *edgelist)
{
    int vertices = getNumberOfVertices(graphString);

    fprintf(stderr,"number_of_vertices: %d\n",vertices);

    if (loadGraph(graphString, vertices, edgelist) == -1)
    {
        return 1;
    }
    return 0;
}

/**
 * a function to transform an edge coloured graph into a vertex coloured simple graph
 */
void generate_expanded_graph(graph *g, int n, int m, int *lab, int*ptn, int vertices, bitset expanded_edgelist){
    EMPTYGRAPH(g, m, n);//
    int edge = 0; // variable to reduce the amount of times the edge-index is computed
    for(int i=0;  i < vertices; i++){ // for every vertex in the graph
            lab[i] = i; //the label of the vertex in the expanded is the label of the vertex in the original graph
            ptn[i] = 1; //all vertices are in the same colour class
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

void compute_possible_edges(struct mygraph *graph_to_label,int* orbits){
    int vertices = graph_to_label->numberOfVertices;
    graph_to_label->possible_edges = EMPTY;
    bitset corrected_orbit = EMPTY;
    for(int edge = graph_to_label->last_edge;  edge < vertices*(vertices-1)/2; edge++){
            if(!contains(corrected_orbit,orbits[vertices+ edge]-vertices) && (edge > graph_to_label->last_edge)){
                add(corrected_orbit,orbits[vertices+ edge]-vertices);
                if(!contains(graph_to_label->edgelist,2*edge)){
                    
                    add(graph_to_label->possible_edges,edge);
                    
                    }
            }
            
    }
}


/**
 * a function to transform a given edgelist of a graph in its canonically labelled form.
 * During this process the graph datatype is also updated with the list of edge-orbits.
 * For each orbit the edge with the lowest index is stored that is higher than the index of the edge that was last added to the graph.
 */
void determine_canonical_labelling(struct mygraph* graph_to_label)
{
    int vertices= graph_to_label->numberOfVertices; //an intermediary variable to store the number of vertices
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
    options.getcanon = TRUE; // enabling nauty to generate a canonical labelling
    options.defaultptn = FALSE; // enabling the option to provide custom colour classes for the vertices
    statsblk stats;// nauty setup
    
    //generate the expanded graph
    generate_expanded_graph(g,n,m,lab,ptn,vertices,graph_to_label->edgelist);
    
    // letting nauty compute the canonical labelling
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);

    // transforming the labelling from the output format of nauty to a map from old label to new label
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
            if (graph_to_label->edgelist & (singleton (2*edge))) {
                add(encoded_canonical,2*mapped_edge);
            }
            if(graph_to_label->edgelist & (singleton (1+2*edge))){
                add(encoded_canonical,2*mapped_edge+1);
            }
            
            edge++;
        }
    }
    
    /**
     * computing the orbits
     */
    compute_possible_edges(graph_to_label,orbits);
   
    //setting the edgelist to the canonically labeled edgelist
    graph_to_label->edgelist= encoded_canonical;
    
    /**
     * freeing temporary variables
     */
    free(lab);
    free(ptn);
    free(mapping);
}


/**
 * A function that generates all non-isomorphic graphs that have start_graph as a sub-graph and
 * have exactly number_to_colour red edges more than start_graph
 * 
 * This function consists of a loop adding a red edge in every edge-orbit that is not coloured yet,
 * the edge in the orbit is always selected as the edge with the smallest (lexicographical) index 
 * larger than the last edge added to the graph.
 */
struct mygraph* generate_final_configuration(struct mygraph* start_graph, int number_to_colour, int* nb_of_graphs)
{

    //initialise a variable to store generated graphs
    struct mygraph* generated_graphs = malloc(sizeof(struct mygraph)*1);
    generated_graphs[0].edgelist = start_graph->edgelist;
    generated_graphs[0].numberOfVertices =start_graph->numberOfVertices;
    start_graph->last_edge = -1;//set the last added edge of the start-graph to -1

    //initialise a variable to store how many graphs were in the previous step
    int nb_of_old_graphs = 1;
    //initialise a variable to store the upper bound on the number of edges that can be coloured in the next step
    int nb_of_empty_edges = 0;

    //determine the canonical labelling of the start-graph and more importantly:
    //determine the list of possible edges to colour
    determine_canonical_labelling(start_graph);

    //the size of this list of edges to colour is the number of empty edges
    nb_of_empty_edges = size(start_graph->possible_edges);
    generated_graphs[0].possible_edges = start_graph->possible_edges;
    generated_graphs[0].last_edge = -1;

    //initialising a variable to store how many new non-isomorphic graphs were generated
    int nb_of_new_graphs = 1;
    //an integer to store how many extra edges were coloured red
    int edges_coloured = 0;
    //a temporary graph variable to store the the new child graph befor it is sorted and inserted in the list
    struct mygraph new_graph;
    

    
    //colour edges red until the required number is reached
    while( number_to_colour > edges_coloured){
        fprintf(stderr,"number of new graphs: %d\n",nb_of_old_graphs);
        
        nb_of_new_graphs = 0;
        
        // the memory allocated for storing the next layer in original labelling with all struct data
        struct mygraph* next_layer_original = malloc(sizeof(struct mygraph)*nb_of_empty_edges);

        //the memory allocated for storing the canoncially labelled next layer
        bitset* next_layer_canonical = malloc(sizeof(bitset)*nb_of_empty_edges);


        nb_of_empty_edges = 0;

        //for each existing graph, colour each of the edges in the set of possible edges
        for (int existing_graph_number = 0; existing_graph_number < nb_of_old_graphs; existing_graph_number ++){
            forEach(edge,generated_graphs[existing_graph_number].possible_edges){
                
            //check if the edge is colourable and not coloured     
            if (!contains(generated_graphs[existing_graph_number].edgelist,2*edge) && contains(generated_graphs[existing_graph_number].edgelist,2*edge+1))
            {
                //make a new graph
                new_graph.numberOfVertices = generated_graphs[existing_graph_number].numberOfVertices;
                
                new_graph.edgelist = generated_graphs[existing_graph_number].edgelist;
                //colour the edge red
                add(new_graph.edgelist,2*edge);
                
                //store the original labelling for later use
                bitset original_edgelist = new_graph.edgelist;

                //set the last added edge to the edge we just coloured
                new_graph.last_edge = edge;

                
                //start by setting possible edges to empty
                new_graph.possible_edges = EMPTY;

                //determine the canonically labeled graph and the set of possible edges
                determine_canonical_labelling(&new_graph);
                
                //add the amount of possible edge to the total amount of possible new graphs for the next iteration
                nb_of_empty_edges += size(new_graph.possible_edges);
                
                //perform a binary search in the list of canonical labellings to check isomorphic copies
                int existing_child = nb_of_new_graphs/2;
                    int begin = 0;
                    int end = nb_of_new_graphs-1;
                    while ((begin <= end) && (next_layer_canonical[existing_child] != new_graph.edgelist))
                    {
                        if (next_layer_canonical[existing_child] < new_graph.edgelist)
                        {
                            begin = existing_child + 1;
                        }
                        else
                        {
                            end = existing_child - 1;
                        }
                        existing_child = (begin + end) / 2;
                    }
                    //we only care about the case where there is no isomorphic copy yet, which is the following
                    if  (existing_child >= nb_of_new_graphs || //if either we have reached the end of the list while searching
                        next_layer_canonical[existing_child]!= new_graph.edgelist) //or we have found where the graph should be but is not
                    {   
                        //next we insert the graph in the correct spot which in the case that 
                        if((nb_of_new_graphs > 0) && //the list is not empty and
                        (next_layer_canonical[existing_child] < new_graph.edgelist))//the search did not end at the end of the list
                        { 
                            //is to the right of where the search ended
                            existing_child +=1;
                        }
                        //copy all graphs in the ordered list right of the new graph one spot to the right
                        for(int j = nb_of_new_graphs; j > existing_child; j--){
                            next_layer_canonical[j] = next_layer_canonical[j-1];
                            
                        }

                        //and add the new graph
                        next_layer_canonical[existing_child] = new_graph.edgelist;

                        //seeing as we never change anything to graphs that have already been generated
                        //the list with the non-canonically labeled graphs is not ordered
                        //the form is used in the next iteration so the notion of last edge remains constant
                        next_layer_original[nb_of_new_graphs].numberOfVertices = new_graph.numberOfVertices;
                        next_layer_original[nb_of_new_graphs].possible_edges = new_graph.possible_edges;
                        next_layer_original[nb_of_new_graphs].edgelist = original_edgelist;
                        next_layer_original[nb_of_new_graphs].last_edge = new_graph.last_edge;
                        nb_of_new_graphs++;
                    }
                }
            }
            
                

            
        }
        //all ways to add this edge were tried, so we move on to the next edge
        edges_coloured ++;

        //free the list of previous graphs
        free(generated_graphs);
        free(next_layer_canonical);
        //the old graphs are now the new graphs
        generated_graphs = next_layer_original;
        nb_of_old_graphs = nb_of_new_graphs;


   }
   *nb_of_graphs = nb_of_old_graphs;
   
   return generated_graphs;
    
}


int main(int argc, char **argv)
{   
    unsigned long vertices = strtol(argv[1], NULL, 10); //the first argument is the number of vertices in the games
    bitset base_graph = singleton(vertices*(vertices-1)/2)-1; //this variable is a bitset with al edges that are allowed to be used
    bitset red_start_graph = EMPTY; //this variable is a bitset with the red subgraph that was present before the players began
    bitset blue_start_graph = EMPTY; //this variable is a bitset with the blue subgraph that was present before the players began
    
    unsigned long bias = 0; //this variable contains if there was a bias between how many edges blue coloured in its turn and how many red coloured

    bool red_player_starts = true; //this variable contains if the red player started the game, if no first player is given red is taken as default
    
    bool red_start_graph_given = false; //this variable stores if the user supplied a value for the red start graph
    bool base_graph_given = false; //this variable stores if the user supplied a value for the base graph
    

    /**
     * this loops over all input arguments to find the optional arguments
     * flags can appear in any order
     * the order of the supplied graphs is fixed: first a base graph, then the red starting graph and then the blue starting graph
     * the first argument is the order, thus we start searching from the second spot
     */
    for(int i =2; i < argc;i++){

        if(argv[i][0] == '-'){//if a flag is detected, check which it is

            if(argv[i][1] == 'h'){//if it was the help flag
                //display the help message
                fprintf(stderr,USAGE);
                exit(-1);
            }

            if(argv[i][1] == 'b'){//if it was the bias flag
                //take the second character as the bias
                bias += strtol(&argv[i][2], NULL, 10);
            }

            if(argv[i][1] == 'p'){//if it was the player flag
                if (strtol(&argv[i][2], NULL, 10)==2){
                    //check if the user wanted player two (blue) to start
                    red_player_starts = false;
                }
            }
        }else{//if an argument was given that is not a flag, assume it is one of the three graph6 string
        if(!base_graph_given){//start with the base graph
            //to read the base graphs in graph6 format they must end with an endofline '\n' character. 
            //this is not parsed correctly, thus these are added in
            char *base_graph_str_complete = malloc(strlen(argv[i]) + 2);
            strcpy(base_graph_str_complete, argv[i]);
            strcat(base_graph_str_complete, "\n");
            
            //then the graph6 string is read
            readGraph(base_graph_str_complete,&base_graph);
            //and the memory for the extended notation freed
            free(base_graph_str_complete);
            base_graph_given = true;

        }else{//if the base graph was set previously
            if(!red_start_graph_given){//check the red starting graph
                //to read the starting graphs in graph6 format they must end with an endofline '\n' character. 
                //this is not parsed correctly, thus these are added in
                char *red_start_graph_str_complete = malloc(strlen(argv[i]) + 2);
                strcpy(red_start_graph_str_complete, argv[i]);
                strcat(red_start_graph_str_complete, "\n");
                
                //then the graph6 string is read
                readGraph(red_start_graph_str_complete,&red_start_graph);
                //and the memory for the extended notation freed
                free(red_start_graph_str_complete);
                red_start_graph_given = true;
    
            }else{//if the red starting graph was already given, assume it is the blue starting graph
                char *blue_start_graph_str_complete = malloc(strlen(argv[i]) + 2);
                strcpy(blue_start_graph_str_complete, argv[i]);
                strcat(blue_start_graph_str_complete, "\n");
    
                readGraph(blue_start_graph_str_complete,&blue_start_graph);
                free(blue_start_graph_str_complete);
    
            }
        }
    }
        

    }
    
    int number_of_graphs;//initialise a variable to store the number of generated graphs

    int number_to_colour = size(base_graph) -size(red_start_graph); //compute the number of red edges in the end configuration 
    //starting from the number of edges that can be coloured and are not yet coloured
    
    if(red_player_starts){
        //if red starts, divide by 2+bias and round up (by adding 1+bias and rounding down)
        number_to_colour = (number_to_colour+1+bias)/(2+bias);
    }else{
        //if blue starts, divide by 2+bias and round down
        number_to_colour = (number_to_colour)/(2+bias);
    }
    fprintf(stderr,"number to colour: %i\n",number_to_colour);
    
    //then create the first graph in the expanded form, using two bits per edge
    bitset start_expanded = EMPTY;
    int edge = 0;//for each edge
    for (int i = 0; i < vertices; i++){
        for (int j =0; j<i; j++){
            if (contains(base_graph,edge)){
                //if it is in the base graph, set its bitvalue to "can be coloured"==01
                removeElement(start_expanded,2*edge);
                add(start_expanded,1+2*edge);

            }
            if (contains(red_start_graph,edge)){
                //if it is in the red starting graph, set its bitvalue to "is coloured red"==11
                //it is assumed here that only edges that are in the base graph are in the starting graphs
                add(start_expanded,2*edge);
                add(start_expanded,1+2*edge);

            }
            if (contains(blue_start_graph,edge)){
                //if it is in the blue starting graph, set its bitvalue to "is coloured blue"==10
                //it is assumed here that only edges that are in the base graph are in the starting graphs
                //it is also assumed here that only edges that are not in the red starting graph are in the blue starting graph
                
                add(start_expanded,2*edge);
                removeElement(start_expanded,1+2*edge);
            }
            edge ++;
        }
    }

    //print this first graph
    fprintf(stderr,"start_Graph_expanded: ");
    for(int j = 0; j < vertices*(vertices-1); j++){
        fprintf(stderr, "%d", (start_expanded & (singleton (j))) ? 1 : 0);
    }
    fprintf(stderr, "\n");

    //add the edgelist to a struct representing a graph
    struct mygraph start_graph;
    start_graph.numberOfVertices = vertices;
    start_graph.edgelist = start_expanded;

    //find the final configurations
    struct mygraph* configurations = generate_final_configuration(&start_graph,number_to_colour,&number_of_graphs);
    //print the number of found configuration so the solver can allocate space
    fprintf(stdout, "%d\n", number_of_graphs);

    //now colour every edge that is not coloured blue to create the final configurations
    for(int conf = 0; conf < number_of_graphs; conf++){//for each graph
        for(int i =0; i < configurations[conf].numberOfVertices;i++){
            for (int j =0; j<i; j++){//for each edge
                if (!contains(configurations[conf].edgelist,2*indexInEdgeList(i,j)) && //if it is not coloured
                contains(configurations[conf].edgelist,2*indexInEdgeList(i,j)+1)){ //and it is an edge that can be coloured
                    //set the state to "is coloured blue"==10
                    removeElement(configurations[conf].edgelist,2*indexInEdgeList(i,j)+1);
                    add(configurations[conf].edgelist,2*indexInEdgeList(i,j));
                }
               
            }
        }  
    }

    //then output all graphs
    for(int conf = 0; conf < number_of_graphs; conf++){

        fprintf(stdout, "%lu\n", configurations[conf].edgelist);
    }
    
    //and free memory used
    free(configurations);

    return 0;
}
