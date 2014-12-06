CC = mpicc
CFLAGS = -Wall

EXEC = average

all: $(EXEC)

average : src/average.c
	$(CC) $(CFLAGS) $^ -o $@ -lm

clean:
	rm -f $(EXEC) *~

.PHONY: clean
