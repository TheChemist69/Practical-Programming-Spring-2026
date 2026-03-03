CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm
TIME = /usr/bin/time --format "$$N %e" --append --output $@

Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : main.o
	$(CXX) $^ -o $@ $(LDLIBS)

main.o : main.cc vector.h matrix.h qr.h lu.h cholesky.h householder.h givens.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

out.times.data : main
	>$@
	for N in $$(seq 100 50 1000); do \
		$(TIME) ./main -size $$N 1>out 2>err ;\
	done
	@cat $@

plot.svg : out.times.data plot_times.gnuplot
	gnuplot plot_times.gnuplot

clean :
	$(RM) main *.o Out.txt out.times.data out err plot.svg

.PHONY : clean
