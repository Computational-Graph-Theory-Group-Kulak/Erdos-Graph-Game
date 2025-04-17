#include <stdio.h>
#include <string.h>
#include "writeGraph6.h"
#include <stdlib.h>
#include "../bitset.h"


char* writeGraph(int n, bitset* adjacencyMatrix) {
    if (n <= 0 || adjacencyMatrix == NULL) {
        return NULL;
    }

    // Calculate the length of the bit vector
    int bitVectorLength = n * (n - 1) / 2;
    int bytesLength = (bitVectorLength + 5) / 6;  // ceil(bitVectorLength / 6)


    // Calculate N(n)
	int headerLength = 1;
	int Nn = n + 63;
    if (n > 62) {
        if (n <= 258047) {
            Nn = 126;
            int tmp = n - 63;
            for (int i = 0; i < 3; ++i) {
                Nn = (Nn << 8) | ((tmp >> (16 - i * 8)) & 0xFF);
            }
			headerLength = 4;
        } else {
            Nn = 126126;
            long long tmp = n - 63;
            for (int i = 0; i < 7; ++i) {
                Nn = (Nn << 8) | ((tmp >> (56 - i * 8)) & 0xFF);
            }
			headerLength = 8;
        }
    }

    // Allocate memory for the graph6 string
    char* graph6String = (char*)malloc((bytesLength + headerLength+1) * sizeof(char)); // Add extra space for header and potential N(n)

    // Write N(n) to the graph6 string
	graph6String[0] = Nn;

    // Write the adjacency matrix bit vector to the graph6 string
    char* bitVector = malloc(bytesLength * sizeof(char));
	for(int i = 0; i < bytesLength; i++) {
		bitVector[i] = 0;
	}
    int bitIndex = 0;
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            int byteIndex = bitIndex / 6;
            int bitOffset = bitIndex % 6;
			if(intersection(adjacencyMatrix[i],singleton(j))){
				bitVector[byteIndex] |= 1 << (5 - bitOffset);
			}

            bitIndex++;
        }
    }

	for(int i = 0; i < bytesLength; i++) {
		bitVector[i] += 63;
	}

	for(int i = 0; i < bytesLength; i++) {
		graph6String[headerLength + i] = bitVector[i];
	}
    free(bitVector);

    // Terminate the graph6 string
    graph6String[headerLength + bytesLength] = '\0';

    return graph6String;
}