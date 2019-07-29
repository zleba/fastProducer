workDir=$PWD 

mkdir -p  nlojetBuild 
wget -qO- https://fastnlo.hepforge.org/code/other/nlojet++-4.1.3-patched.tar.gz    | tar zxv 
mv nlojet++-4.1.3  nlojet

cd nlojet
./configure --prefix=$workDir/nlojetBuild
make -j`nproc`
make install
cd ..
