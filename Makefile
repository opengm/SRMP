DIRS := src/examples src/FactorTypes src/Algs

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CFLAGS := -O3 -DNDEBUG
#CFLAGS := -g
CXX ?= c++
LIBS := 
INCLUDES := -I./include
LIBDIR := 

# Add librt if the target platform is not Darwin (OS X)
ifneq ($(shell uname -s),Darwin)
    LIBS += -lrt
endif
 
all: srmp

srmp: ${OBJS}
	$(CXX) $(CFLAGS) ${LIBDIR} -o $@ ${OBJS} ${LIBS}

.cpp.o:
	$(CXX) $(CFLAGS) ${INCLUDES} $< -c -o $@

clean:
	rm -f ${OBJS} srmp
