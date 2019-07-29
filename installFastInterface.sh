workDir=$PWD 
mkdir -p  fastnloInterfaceBuild  

wget -qO- https://fastnlo.hepforge.org/code/v23/fastnlo_interface_nlojet-2.3.1pre-2411.tar.gz  | tar zxv 
mv fastnlo_interface_nlojet-2.3.1pre-2411 fastnloInterface

cd fastnloInterface
./configure --prefix=$workDir/fastnloInterfaceBuild --with-nlojet=$workDir/nlojetBuild --with-fnlotoolkit=$workDir/fastnloBuild
make -j`nproc`
make install
cd ..
