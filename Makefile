all: ssss-split

ssss-split: ssss-split.o

ssss-split.o: vendor/ssss.c
	$(CC) -c -O2 $< -o $@

ssss-split: ssss-split.o
	$(CC) -O2 $< -o $@ -lgmp

clean:
	$(RM) ssss-split ssss-split.o
