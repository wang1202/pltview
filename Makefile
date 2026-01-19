# Makefile for pltview (C version)

CC = gcc
CFLAGS = -O3 -Wall -march=native
LDFLAGS = -lX11 -lXt -lXaw -lXmu -lm

# macOS specific
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CFLAGS += -I/opt/X11/include
    LDFLAGS += -L/opt/X11/lib
endif

TARGET = pltview_c
SRC = pltview.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

.PHONY: all clean install
