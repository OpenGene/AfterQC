TARGET = editdistance/libed.so
CC = g++
${TARGET}:${OBJ}
	$(CC) editdistance/_editdistance.cpp -fPIC -shared -O3 -o editdistance/libed.so