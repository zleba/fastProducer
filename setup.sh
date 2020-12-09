if [ -d "/cvmfs/sft.cern.ch/lcg" ]; then
    source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-centos7-gcc7-opt/setup.sh
    export LHAPDF_DATA_PATH=$PWD/dpdfLib/DPDFsets:/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:`lhapdf-config --datadir`
else
    echo "The cvmfs file system not mounted"
fi

export PlH_DIR=$PWD/PlottingHelper
export PROJECT_DIR=$PWD

LD_LIBRARY_PATH=$PWD/nlojetBuild/lib:$LD_LIBRARY_PATH
PATH=$PWD/nlojetBuild/bin:$PATH
PATH=$PWD/fastnloBuild/bin/:$PATH
