CC=mpicc
EXEC= main.out

all: $(EXEC) 

%.out: %.o
	$(CC) $< -o $@ -lm

%.o: %.c
	$(CC) -c $< -lm

clean: 
	rm -f *.out sol* time* speedup