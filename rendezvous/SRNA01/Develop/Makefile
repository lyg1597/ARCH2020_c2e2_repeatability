test: annotation.o point.o rep-point.o linear-set.o reach-tube.o checker.o simulator.o generalized-star.o gstar-reach-tube.o verifier.o model.o main.o
	g++ -march=native -std=c++11 annotation.o point.o rep-point.o linear-set.o reach-tube.o checker.o simulator.o generalized-star.o gstar-reach-tube.o verifier.o model.o main.o -g -I /usr/include/python3.5m -lpython3.5m -lppl -lglpk -o test -lgmp -ldl

annotation.o: annotation.cpp annotation.hpp
	g++ -march=native -std=c++11 -c annotation.cpp -g

point.o: point.cpp point.hpp
	g++ -march=native -std=c++11 -c point.cpp -g

rep-point.o: rep-point.cpp rep-point.hpp
	g++ -march=native -std=c++11 -c rep-point.cpp -g

linear-set.o: linear-set.cpp linear-set.hpp
	g++ -march=native -std=c++11 -c linear-set.cpp -g

reach-tube.o: reach-tube.cpp reach-tube.hpp
	g++ -march=native -std=c++11 -I ../eigen-eigen-323c052e1731 -c reach-tube.cpp -g 

checker.o: checker.cpp checker.hpp
	g++ -march=native -std=c++11 -c checker.cpp -g

simulator.o: simulator.cpp simulator.hpp
	g++ -march=native -std=c++11 -I /usr/include/python3.5m -c simulator.cpp -g 

generalized-star.o: generalized-star.cpp generalized-star.hpp
	g++ -march=native -std=c++11 -I ../eigen-eigen-323c052e1731 -c generalized-star.cpp -g 

gstar-reach-tube.o: gstar-reach-tube.cpp gstar-reach-tube.hpp
	g++ -march=native -std=c++11 -I ../eigen-eigen-323c052e1731 -c gstar-reach-tube.cpp -g 

verifier.o: verifier.cpp verifier.hpp
	g++ -march=native -std=c++11 -I ../eigen-eigen-323c052e1731 -I /usr/include/python3.5m -c verifier.cpp -g

model.o: model.cpp model.hpp
	g++ -march=native -std=c++11 -I ../eigen-eigen-323c052e1731 -I /usr/include/python3.5m -c model.cpp -g

main.o: main.cpp
	g++ -march=native -std=c++11 -I ../eigen-eigen-323c052e1731 -c main.cpp -g

clean: 
	rm *.o test