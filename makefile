compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.

solver-64bit: algorithm-erdos-multi-nauty.c readGraph/readGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_64_BIT -o ErdosGame-multi-nauty algorithm-erdos-multi-nauty.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c -DUSE_TLS $(flags) -O3

solver-64bit-gprof: algorithm-erdos-multi-nauty.c readGraph/readGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_64_BIT -o ErdosGame-multi-nauty algorithm-erdos-multi-nauty.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c -DUSE_TLS $(flags) -g -pg -fsanitize=address

generator-64bit-gprof: gen-colex.c readGraph/readGraph6.c writeGraph/writeGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_64_BIT -o gen-colex gen-colex.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c writeGraph/writeGraph6.c -DUSE_TLS $(flags) -g -pg -fsanitize=address

generator-64bit: gen-colex.c readGraph/readGraph6.c writeGraph/writeGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_64_BIT -o gen-colex gen-colex.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c writeGraph/writeGraph6.c -DUSE_TLS $(flags) -O3

solver-128bit: algorithm-erdos-multi-nauty.c readGraph/readGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_128_BIT -o ErdosGame-multi-nauty algorithm-erdos-multi-nauty.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c -DUSE_TLS $(flags) -O3

solver-128bit-gprof: algorithm-erdos-multi-nauty.c readGraph/readGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_128_BIT -o ErdosGame-multi-nauty algorithm-erdos-multi-nauty.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c -DUSE_TLS $(flags) -g -pg -fsanitize=address

generator-128bit-gprof: gen-colex.c readGraph/readGraph6.c writeGraph/writeGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_128_BIT -o gen-colex gen-colex.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c writeGraph/writeGraph6.c -DUSE_TLS $(flags) -g -pg -fsanitize=address

generator-128bit: gen-colex.c readGraph/readGraph6.c writeGraph/writeGraph6.c bitset.h nauty2_8_6/nauty.h
	$(compiler) -lpthread -DUSE_128_BIT -o gen-colex gen-colex.c  nauty2_8_6/nautyL1.a -DMAXN=WORDSIZE readGraph/readGraph6.c writeGraph/writeGraph6.c -DUSE_TLS $(flags) -O3

all-64bit: generator-64bit solver-64bit

all-128bit: generator-128bit solver-128bit

.PHONY: clean
clean:
	rm -f 

