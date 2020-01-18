# Makefile - basic makefile rules for testing the header

.PHONY: default clean

default: test.out

clean:
	rm -rf $(wildcard test.out)

%.out: %.c cd_vector.h
	$(CC) $< -lm -o $@

