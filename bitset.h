#ifndef BITSETCHOOSER
#define BITSETCHOOSER

#ifdef USE_8_BIT
	#include "bitset8Vertices.h"
	#define MAXVERTICES 8

#elif USE_64_BIT
	#include "bitset64Vertices.h"
	#define MAXVERTICES 64

#elif defined(USE_128_BIT)
	#include "bitset128Vertices.h"
	#define MAXVERTICES 128

#elif defined(USE_128_BIT_ARRAY)
	#include "bitset128VerticesArray.h"
	#define MAXVERTICES 128

#endif

#endif
