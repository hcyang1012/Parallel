all:
	gcc -o 2mm 2mm.c -g
	gcc -o 2mm_vector 2mm.c -DCONFIG_SSE2 -DCONFIG_VECTOR -msse2 -g -O2

test:
	./2mm 1 10 10 10 10 2> test.out
	./2mm_vector 1 10 10 10 10 2> test2.out
	diff test.out test2.out

perf:
	time ./2mm 0 1200 1200 1200 1200 
	time ./2mm_vector 0 1200 1200 1200 1200

clean:
	rm 2mm 2mm_vector
