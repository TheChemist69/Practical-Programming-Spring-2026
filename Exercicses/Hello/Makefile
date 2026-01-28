Out.txt : main            
	./main > Out.txt    


main : main.o				 
	c++ -o main main.o

main.o: main.cc                 
	c++ -c main.cc -o main.o


clean:
	rm -f Out.txt main

test:
	echo CXX
	echo $(CXX)
	