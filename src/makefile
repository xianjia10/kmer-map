object = main.o kmerfind.o generatemap.o
target = kmermap
kmermap: $(object)
	g++ -o $(target) $(object) -lpthread
	mv $(target) ../
.PHONY : clean
clean :
	rm $(objects)
