#CXX=g++-4.7
CXX=clang++
CXXFLAGS=-std=c++11
CXXINCLUDE=-I../
CXXLDFLAGS=-lm -lpthread -ldl -lcrypt

OUT=a.out
CGOUT=cg.out

$(OUT) : main.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) $(CXXLDFLAGS) -o $(OUT) $< 

run: $(OUT)
	./$(OUT)

codingame : $(OUT)
	$(shell cat ../algo/detail/Node.hpp ../algo/detail/GraphSearchHelpers.hpp ../algo/Graph.hpp ../algo/GraphSearch.hpp ./main.cpp > codingame.cpp)
	$(shell sed -i "s/#include \".*//g" codingame.cpp)
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) $(CXXLDFLAGS) -o $(CGOUT) codingame.cpp
	xclip -selection clipboard -in codingame.cpp

clean:
	rm -f $(OUT)
	rm -f $(CGOUT)
	rm -f codingame.cpp
