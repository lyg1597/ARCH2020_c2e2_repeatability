g++ -g -w -O2 -std=c++11 simulator.cpp -o simu
g++ -g -fPIC -shared hybridSimGI.cpp -o libhybridsim.so -lppl -lgmp
g++ -g -fPIC -shared bloatedSimGI.cpp -o libbloatedsim.so -lppl -lgmp