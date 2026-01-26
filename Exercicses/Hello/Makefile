Out.txt : hello             # Out.txt depends on hello.exe
	./hello > Out.txt        # run hello.exe, send output to Out.txt
hello: hello.cc                 # hello depends on hello.cc
	c++ -o hello hello.cc  # compile-and-link hello.cc, save executable in hello
clean:
	rm -f Out.txt hello