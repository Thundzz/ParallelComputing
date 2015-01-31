CC=mpicc
EXEC= main.out

all: $(EXEC) 

%.out: %.o
	$(CC) -W -Wall $< -o $@ -lm

%.o: %.c
	$(CC) -W -Wall -c $< -lm

clean: 
	rm -f *.out sol* time* speedup