CC = mpicc
CFLAGS = -Wall

EXEC = average constants sparse starter

all: $(EXEC)

average : src/gfx.c src/average.c
	$(CC) $(CFLAGS) $^ -o $@ -lm -lX11

starter : src/starter.c
	$(CC) $(CFLAGS) $^ -o $@ -lm

constants : src/gfx.c  src/constants.c
	$(CC) $(CFLAGS) $^ -o $@ -lm -lX11

gfx : src/gfx.c src/example.c
	gcc $(CFLAGS) $^ -o $@ -lm -lX11

sparse : src/gfx.c  src/sparse.c
	$(CC) $(CFLAGS) $^ -o $@ -lm -lX11

clean:
	rm -f $(EXEC) *~
	rm -f src/*~

.PHONY: clean
