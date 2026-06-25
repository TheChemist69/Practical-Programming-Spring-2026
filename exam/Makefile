CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lm

OBJS = main.o cintegrator.o bessel.o

PLOTS = plot_bessel.svg plot_bessel_convergence.svg \
        plot_real_trig.svg plot_real_rectangle.svg \
        plot_spectral.svg

all: Out.txt $(PLOTS)

# Build and run; the program writes Out.txt itself (and echoes to the console).
Out.txt: main
	./main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)

# All data files are produced together by the single run.
DATA = bessel_contour.data bessel_polygon_convergence.data \
       real_trig.data real_rectangle.data \
       spectral_trapezoid.data spectral_adaptive.data
$(DATA): Out.txt

main.o: main.cc bessel.h cintegrator.h
	$(CXX) $(CXXFLAGS) -c main.cc -o $@

cintegrator.o: cintegrator.cc cintegrator.h
	$(CXX) $(CXXFLAGS) -c cintegrator.cc -o $@

bessel.o: bessel.cc bessel.h cintegrator.h
	$(CXX) $(CXXFLAGS) -c bessel.cc -o $@

# --- plots (skipped gracefully if gnuplot is absent) ---
plot_%.svg: plot_%.gnuplot $(DATA)
	@if command -v gnuplot >/dev/null 2>&1; then gnuplot $<; \
	else echo "gnuplot not found, skipping $@"; fi

clean:
	$(RM) main $(OBJS) Out.txt *.data *.svg

.PHONY: all clean
