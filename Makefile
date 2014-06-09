CC = gcc
CXX = g++

CXXFLAGS = -c -g -Wall -Wextra -pedantic 
CCFLAGS = -c
LDFLAGS = -lpthread

CPP_SRC = \
	main.cpp \
	pdbhelper.cpp \
	pdbatom.cpp

CPP_SRC_OBJS = $(CPP_SRC:.cpp=.o)

CC_SRC = file_read_write.c
CC_SRC_OBJS = $(CC_SRC:.c=.o)

EXECUTABLE = main

$(EXECUTABLE): $(CPP_SRC_OBJS) $(CC_SRC_OBJS)
	$(CXX) $(CPP_SRC_OBJS) $(CC_SRC_OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

.c.o:
	$(CC) $(CCFLAGS) $< -o $@

all: $(SOURCES) $(EXECUTABLE)
	
clean:
	rm -rf $(CPP_SRC_OBJS) $(CC_SRC_OBJS) $(EXECUTABLE)
