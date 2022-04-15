object = main.o kmerfind.o gen_map.o
target = kmermap
kmfilter: $(object)
	g++ -o $(target) $(object) -lpthread
	mv $(target) ../
.PHONY : clean
clean :
	rm $(objects)
