# Compiler and flags
CC = gcc

CFLAGS = -O2 -g -march=x86-64 -funroll-loops -ffast-math -Wall -Wextra -pthread -I.

LDFLAGS = -lm -pthread

# Object files
OBJ = virumble.o hash.o mem.o merge.o aux.o

# Default target: build the executable
.PHONY: all
all: virumble

# Rule to compile any .c file into .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link final binary with all object files
virumble: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

# Run the program
run: virumble
	./virumble

# Profile with gprof (still builds all objects with profiling flags)
profile: 
	$(MAKE) CFLAGS="-pg -g -O0 -march=native -Wall -Wextra -pthread -I." \
	        LDFLAGS="-pg -lm -pthread" \
	        virumble
	./virumble
	gprof virumble gmon.out > profile.txt
	@echo "Profile report saved to profile.txt"

# Clean build artifacts
.PHONY: clean
clean:
	rm -f *.o virumble dataset_modifier gmon.out profile.txt