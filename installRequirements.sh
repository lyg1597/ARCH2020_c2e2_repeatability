#!/bin/bash

# APT_IN="apt-get install -y"
# APT_IN_NR="apt-get install -y --no-install-recommends"

# if [ $UID != 0 ]
# then
#   echo "Please run this script with sudo: sudo $0"
#   exit 1
# fi

echo "Dependencies"
apt-get update
apt-get upgrade
apt-get install -y python3-dev
apt-get install -y build-essential

apt-get remove -y glpk
apt-get remove -y libglpk-dev
apt-get remove -y libglpk0
wget -O glpk-4.55.tar.gz http://ftp.gnu.org/gnu/glpk/glpk-4.55.tar.gz
tar -xvhf glpk-4.55.tar.gz
cd glpk-4.55
./configure
make
make install
cd ..

apt-get install -y libglpk-dev
apt-get install -y libglpk0
apt-get install -y bison
wget -O bison-3.0.2.tar.gz http://ftp.gnu.org/gnu/bison/bison-3.0.2.tar.gz
tar -xhvf bison-3.0.2.tar.gz
cd bison-3.0.2
./configure
make
make install
cd ..

wget ftp://ftp.cs.unipr.it/pub/ppl/releases/1.2/ppl-1.2.tar.gz
tar -xhvf ppl-1.2.tar.gz
cd ppl-1.2
./configure
make
make install
cd ..

apt-get install -y libboost-dev
apt-get install -y libboost-python-dev



# #python 
# apt-get install -y python3-pip
# apt-get remove -y python3-numpy
# python3 -m pip install --upgrade pip
# # pip3 install sympy
# # pip3 install numpy
# # apt-get install -y python3-tk
# # pip3 install scipy
# # apt-get install -y python3-ply
# apt-get install -y libtiff4-dev libjpeg8-dev zlib1g-dev libfreetype6-dev liblcms2-dev libwebp-dev tcl8.5-dev tk8.5-dev python-tk
# # pip3 install Pillow
# apt-get install -y git

# # apt-get install python3-pil.imagetk


# git clone https://github.com/oblalex/gnuplot.py-py3k.git
# cd gnuplot.py-py3k/
# python3 setup.py install
# cd ..
# apt-get install -y python-gnuplot

# # Bokeh plotting requirements

# # pip3 install bokeh
# # pip3 install selenium
# # apt install -y nodejs-legacy
# # apt install -y npm
# # npm install -g phantomjs-prebuilt

# # Eigen library for backends

# # wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
# # tar xjf 3.3.7.tar.bz2

# #Removes compressed files

# rm *.tar.gz
# rm *.zip
