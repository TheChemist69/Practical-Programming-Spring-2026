CXXFLAGS = -Wall -Wextra -std=c++23 -O
LDLIBS = -lstdc++ -lm -lpthread
N = 1000000000
TIME = /usr/bin/time --format "$$nthreads %e %U" --append --output $@

main : main.o harm.o
	$(CXX) $^ -o $@ $(LDLIBS)

main.o : main.cc harm.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

harm.o : harm.cc harm.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Out.txt : main
	./main > Out.txt
	@cat Out.txt

out.times : main
	>$@
	>log.main
	for nthreads in 1 2 3 4 5 6 7 8; do \
		$(TIME) ./main -nthreads $$nthreads -nterms $(N) >>log.main ;\
	done
	@cat $@

plot.svg : out.times plot_times.gnuplot
	gnuplot plot_times.gnuplot

clean :
	$(RM) main *.o Out.txt out.times log.main plot.svg

.PHONY : clean
