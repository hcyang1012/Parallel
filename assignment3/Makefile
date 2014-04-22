all: 2mm.exe 2mm_tile.exe

2mm.exe: 2mm.c
	gcc 2mm.c -o 2mm.exe

2mm_tile.exe: 2mm.c
	gcc 2mm.c -o 2mm_tile.exe -DCONFIG_TILE -DN=25

test:
	./2mm.exe 1 50 50 50 50 2> test.out -g
	./2mm_tile.exe 1 50 50 50 50 2> test2.out
	diff test.out test2.out

perf:
	time ./2mm.exe 0 500 500 500 500
	time ./2mm_tile.exe 0 500 500 500 500

clean:
	rm 2mm.exe 2mm_tile.exe test.out test2.out
