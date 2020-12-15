wget -qO-  http://ekpwww.etp.kit.edu/~rabbertz/fastNLO_CMS/fastNLO-with-NLOJet++/fastnlo_toolkit-2.3.1-2771.tar.gz | tar zxv
mkdir -p fastnloBuild
mv fastnlo_toolkit-2.3.1-2771 fastnlo
pwd=$PWD
cd fastnlo
./configure --prefix=$pwd/fastnloBuild  --disable-doxygen-doc  --with-lhapdf=$(dirname $(dirname `which lhapdf-config`)) CXXFLAGS=" -O3 -g -fopenmp "
make -j`nproc`
make install
