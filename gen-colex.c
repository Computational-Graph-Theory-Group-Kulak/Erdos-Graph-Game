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

#define indexInEdgeList(vertex_1, vertex_2) (((vertex_1) < (vertex_2)) ? ((vertex_2) * ((vertex_2) - 1) / 2 + (vertex_1)) : ((vertex_1) * ((vertex_1) - 1) / 2 + (vertex_2)))

struct mygraph
{
    int numberOfVertices;
    bitset edgelist;
    bitset orbits;
    int last_edge;

};

struct thread_data
{
    struct mygraph *generated_graphs;
    int nb_of_old_graphs;
    struct mygraph *new_layer;
    int nb_of_new_graphs;
};

struct mergesort_thread_data
{
    struct mygraph *children_1;
    int childrennumber_1;
    struct mygraph *children_2;
    int childrennumber_2;
    struct mygraph *children;
    int *childrennumber;
};


int readGraph(const char *graphString, struct mygraph *g)
{
    g->numberOfVertices = getNumberOfVertices(graphString);

    fprintf(stderr,"number_of_vertices: %d\n",g->numberOfVertices);

    if (loadGraph(graphString, g->numberOfVertices, &g->edgelist) == -1)
    {
        return 1;
    }
    return 0;
}

void determine_canonical_labeling(struct mygraph* graph_to_label)
{
    int vertices= graph_to_label->numberOfVertices;
    int* mapping = malloc(sizeof(int)*vertices);
    int n = vertices+(vertices*(vertices-1)/2)+4;
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
    graph g[n*m];
    graph cg[n*m];
    int * lab = malloc(n * sizeof(int));
    int * ptn = malloc(n * sizeof(int));
    int orbits[n];
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    EMPTYGRAPH(g, m, n);
    int edge = 0;
    for(int i=0;  i < vertices; i++){
            lab[i] = i;
            ptn[i] = 1;
        for(int j = 0; j < i; j++){
            lab[vertices+edge] = vertices+edge;
            ptn[vertices+edge] = 1;
            ADDONEEDGE(g, i, vertices+edge, m);
            ADDONEEDGE(g, j, vertices+edge, m);
            if(contains(graph_to_label->edgelist,2*edge)&&contains(graph_to_label->edgelist,1+2*edge)){
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2)+1, m);
            }else if(contains(graph_to_label->edgelist,singleton (2*edge))){
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2)+2, m);
            }else if (!contains(graph_to_label->edgelist ,1+2*edge)){
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2), m);
            } 
            else{
                ADDONEEDGE(g, vertices+edge, vertices+(vertices*(vertices-1)/2)+3, m);
            }
            edge++;
        }
    }


    lab[vertices+(vertices*(vertices-1)/2)] = vertices+(vertices*(vertices-1)/2);
    lab[vertices+(vertices*(vertices-1)/2)+1] = vertices+(vertices*(vertices-1)/2)+1;
    lab[vertices+(vertices*(vertices-1)/2)+2] = vertices+(vertices*(vertices-1)/2)+2;
    lab[vertices+(vertices*(vertices-1)/2)+3] = vertices+(vertices*(vertices-1)/2)+3;
    ptn[vertices+(vertices*(vertices-1)/2)] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)+1] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)+2] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)+3] = 0;
    ptn[vertices+(vertices*(vertices-1)/2)-1] = 0;
    ptn[vertices-1] = 0;
    
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);

    for (int i = 0; i < vertices; i++) {
            mapping[lab[i]] = i;
    }

    bitset encoded_canonical = EMPTY;
    edge = 0;
    for (int i = 0; i < vertices; i++){
        for (int j = 0; j < i; j++){
            int mapped_edge = indexInEdgeList(mapping[i],mapping[j]);
            if (graph_to_label->edgelist & (singleton (2*edge))) {
                encoded_canonical |= (singleton (2*mapped_edge));
            }
            if(graph_to_label->edgelist & (singleton (1+2*edge))){
                encoded_canonical |= (singleton (2*mapped_edge+1));
            }
            
            edge++;
        }
    }
    graph_to_label->orbits = EMPTY;
    bitset corrected_orbit = EMPTY;
    edge = 0;
    for(int i=0;  i < vertices; i++){
        for(int j = 0; j < i; j++){
            //fprintf(stderr,"vertex map %d vs edge maps %d\n",edge,edge_mapping[indexInEdgeList(i,j)]);
            //fprintf(stderr,"edge: %d,%d mapped edge: %d, orbits with mapping of vertices: %d, orbits with mapping of edge: %d\n",i,j,edge,orbits[vertices+edge]-vertices,orbits[vertices+indexInEdgeList(i,j)]-vertices);
            if(!contains(corrected_orbit,orbits[vertices+ edge]-vertices) && (edge > graph_to_label->last_edge)){
                add(corrected_orbit,orbits[vertices+ edge]-vertices);
                if(!contains(graph_to_label->edgelist,2*edge)){
                    
                    add(graph_to_label->orbits,edge);
                    
                    }
            }

            edge++;
            
        }
    }
   
    graph_to_label->edgelist= encoded_canonical;
    
    free(lab);
    free(ptn);
    free(mapping);
}

struct mygraph* generate_final_configuration(struct mygraph* start_graph, int number_to_color, int* nb_of_graphs, int max_threads)
{

 
    struct mygraph* generated_graphs = malloc(sizeof(struct mygraph)*1);
    int nb_of_old_graphs = 1;
    int nb_of_empty_edges = 0;
    generated_graphs[0].edgelist = start_graph->edgelist;
    generated_graphs[0].numberOfVertices =start_graph->numberOfVertices;
    start_graph->last_edge = -1;
    


    determine_canonical_labeling(start_graph);

    nb_of_empty_edges = size(start_graph->orbits);
    generated_graphs[0].orbits = start_graph->orbits;
    generated_graphs[0].last_edge = -1;
    int nb_of_new_graphs = 1;
    int edges_colored = 0;
    struct mygraph new_graph;
    

    

    while( number_to_color > edges_colored){
        fprintf(stderr,"number of new graphs: %d\n",nb_of_old_graphs);
        
        nb_of_new_graphs = 0;
        struct mygraph* next_layer_original = malloc(sizeof(struct mygraph)*nb_of_empty_edges);
        bitset* next_layer_canonical = malloc(sizeof(bitset)*nb_of_empty_edges);
        nb_of_empty_edges = 0;
        for (int existing_graph_number = 0; existing_graph_number < nb_of_old_graphs; existing_graph_number ++){
            forEach(edge,generated_graphs[existing_graph_number].orbits){
                
                
            if (!contains(generated_graphs[existing_graph_number].edgelist,2*edge) && contains(generated_graphs[existing_graph_number].edgelist,2*edge+1))
            {
                new_graph.numberOfVertices = generated_graphs[existing_graph_number].numberOfVertices;
                
                new_graph.edgelist = generated_graphs[existing_graph_number].edgelist;
                add(new_graph.edgelist,2*edge);

                bitset original_edgelist = new_graph.edgelist;
                new_graph.last_edge = edge;

                

                new_graph.orbits = EMPTY;
                determine_canonical_labeling(&new_graph);
                
                nb_of_empty_edges += size(new_graph.orbits);
                
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
        
                    if  (existing_child < nb_of_new_graphs && next_layer_canonical[existing_child]== new_graph.edgelist)
                        {

                        }
            
                    else
                    {   
                        if((nb_of_new_graphs > 0) && (next_layer_canonical[existing_child] < new_graph.edgelist)){
                            existing_child +=1;
                        }
                        for(int j = nb_of_new_graphs; j > existing_child; j--){
                            next_layer_canonical[j] = next_layer_canonical[j-1];
                            
                        }
                        next_layer_canonical[existing_child] = new_graph.edgelist;
                        next_layer_original[nb_of_new_graphs].numberOfVertices = new_graph.numberOfVertices;
                        next_layer_original[nb_of_new_graphs].orbits = new_graph.orbits;
                        next_layer_original[nb_of_new_graphs].edgelist = original_edgelist;
                        next_layer_original[nb_of_new_graphs].last_edge = new_graph.last_edge;
                        nb_of_new_graphs++;
                    }
                }
            }
            
                

            
        }
    
        edges_colored ++; 
        free(generated_graphs);
        free(next_layer_canonical);
        generated_graphs = next_layer_original;
        nb_of_old_graphs = nb_of_new_graphs;


   }
   *nb_of_graphs = nb_of_old_graphs;
   
   return generated_graphs;
    
}


int main(int argc, char **argv)
{
    char* base_graph_str = argv[1];
    char* start_graph_str = argv[2];
    char* second_start_graph_str = argv[3];
    unsigned long bias = strtol(argv[4], NULL, 10);
    unsigned long threads = strtol(argv[5], NULL, 10);
    
    char *base_graph_str_complete = malloc(strlen(base_graph_str) + 2);
    strcpy(base_graph_str_complete, base_graph_str);
    strcat(base_graph_str_complete, "\n");
    struct mygraph base_graph;
    base_graph.edgelist = EMPTY;
    readGraph(base_graph_str_complete, &base_graph);
    
    char *start_graph_str_complete = malloc(strlen(start_graph_str) + 2);
    strcpy(start_graph_str_complete, start_graph_str);
    strcat(start_graph_str_complete, "\n");
    struct mygraph start_graph;
    start_graph.edgelist = EMPTY;
    readGraph(start_graph_str_complete,&start_graph);

    char *second_start_graph_str_complete = malloc(strlen(second_start_graph_str) + 2);
    strcpy(second_start_graph_str_complete, second_start_graph_str);
    strcat(second_start_graph_str_complete, "\n");
    struct mygraph second_start_graph;
    second_start_graph.edgelist = EMPTY;
    readGraph(second_start_graph_str_complete,&second_start_graph);

    bool first_player = true;
    if (argc > 6){
        if (strtol(argv[6], NULL, 10)>1){
            first_player = false;
        }
        
    }else{

        if(size(start_graph.edgelist)>0){
            if(size(second_start_graph.edgelist) != 0){
                fprintf(stderr,"missing first player\n");
                exit(-1);
            }
            first_player = false;
                
        }
    }
    

    int number_of_graphs;

    int number_to_color = size(base_graph.edgelist) -size(start_graph.edgelist);
    if(first_player){
        number_to_color = (number_to_color+1+bias)/(2+bias);
    }else{
        number_to_color = (number_to_color)/(2+bias);
    }
    
    //number_to_color = 5;
    //fprintf(stderr,"number to color:%d\n",number_to_color);
    bitset start_expanded = EMPTY;
    int edge = 0;
    for (int i = 0; i < base_graph.numberOfVertices; i++){
        for (int j =0; j<i; j++){
            if (base_graph.edgelist & singleton(edge)){
                start_expanded |= (singleton (1+2*edge));

            }
            if (start_graph.edgelist & singleton(edge)){
                start_expanded |= (singleton (2*edge));

            }
            if (second_start_graph.edgelist & singleton(edge)){
                start_expanded |= (singleton (2*edge));
                start_expanded &= ~(singleton (1+2*edge));
            }
            edge ++;
        }
    }
    start_graph.edgelist = start_expanded;
    fprintf(stderr,"start_Graph_expanded: ");
    for(int j = 0; j < start_graph.numberOfVertices*(start_graph.numberOfVertices-1); j++){
        fprintf(stderr, "%d", (start_graph.edgelist & (singleton (j))) ? 1 : 0);
    }
    fprintf(stderr, "\n");
    struct mygraph* configurations = generate_final_configuration(&start_graph,number_to_color,&number_of_graphs, threads);
    fprintf(stdout, "%d\n", number_of_graphs);
    bitset* configuration_adjacency = malloc(sizeof(bitset)*configurations[0].numberOfVertices);
    for(int conf = 0; conf < number_of_graphs; conf++){
        
        for(int i =0; i < configurations[conf].numberOfVertices;i++){
            configuration_adjacency[i] = EMPTY;
        }
        for(int i =0; i < configurations[conf].numberOfVertices;i++){
            
            for (int j =0; j<i; j++){
                if(configurations[conf].edgelist & singleton (2*indexInEdgeList(i,j))){
                    configuration_adjacency[i] |= singleton(j);
                    configuration_adjacency[j] |= singleton(i);
                    //fprintf(stderr, "%d,%d ", i,j);
                }else if (configurations[conf].edgelist & singleton (2*indexInEdgeList(i,j)+1)){
                    configurations[conf].edgelist &= ~singleton (2*indexInEdgeList(i,j)+1);
                    configurations[conf].edgelist |= singleton (2*indexInEdgeList(i,j));
                }
               
            }
        }  
    }
    free(configuration_adjacency);
    for(int conf = 0; conf < number_of_graphs; conf++){

        fprintf(stdout, "%lu\n", configurations[conf].edgelist);
    }
    
    free(configurations);

    return 0;
}