all: 2mm.exe 2mm_tile.exe

2mm.exe: 2mm.c
	gcc 2mm.c -o 2mm.exe

2mm_tile.exe: 2mm.c
	gcc 2mm.c -o 2mm_tile.exe -DCONFIG_TILE

test:
	./2mm.exe 1 10 10 10 10 2> test.out
	./2mm_tile.exe 1 10 10 10 2> test2.out
	diff test.out test2.out

perf:
	time ./2mm.exe 0 100 100 100 100
	time ./2mm_tile.exe 0 100 100 100 100

clean:
	rm 2mm.exe 2mm_tile.exe test.out test2.out
