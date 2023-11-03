#!/bin/bash

wget https://raw.githubusercontent.com/Jonathan-123/links2/master/links2.c
wget https://raw.githubusercontent.com/Jonathan-123/links2/master/links2.h

gcc -c -Wall -Werror -fPIC links2.c -o blah.o
gcc -shared -o links2.so blah.o
rm blah.o

mv links2.h ./usr/include/links2.h
mv links2.so ./usr/lib/liblinks2.so
