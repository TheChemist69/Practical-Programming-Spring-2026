CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++14
TARGET = main

all: run



OBJS = main.o erf_func.o sgamma_func.o lngamma_func.o
lngamma_func.o: lngamma_func.cc lngamma_func.h
	$(CXX) $(CXXFLAGS) -c lngamma_func.cc
sgamma_func.o: sgamma_func.cc sgamma_func.h
	$(CXX) $(CXXFLAGS) -c sgamma_func.cc

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) -lm

main.o: main.cc erf_func.h
	$(CXX) $(CXXFLAGS) -c main.cc

erf_func.o: erf_func.cc erf_func.h
	$(CXX) $(CXXFLAGS) -c erf_func.cc

run: $(TARGET)
	./$(TARGET)
	@echo "Run complete. Data files: erf_curve.dat, erf_points.dat, gamma_curve.dat, gamma_points.dat, lngamma_curve.dat, lngamma_points.dat"

clean:
	rm -f $(TARGET) main.o erf_func.o sgamma_func.o lngamma_func.o erf_curve.dat erf_points.dat gamma_curve.dat gamma_points.dat lngamma_curve.dat lngamma_points.dat