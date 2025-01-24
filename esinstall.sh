#~ git clone -b R_ensemble_merge_python https://github.com/jonaslandsgesell/espresso.git
git clone https://github.com/espressomd/espresso.git

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
#Now close and re-open your terminal window for the changes to take effect. Or say

export PATH="$HOME/miniconda3/bin:$PATH"
export PATH="$HOME/miniconda2/bin:$PATH"

conda install cmake numpy boost sphinx ipython cython hdf5  numpydoc
#~ conda install cython=0.25.2
conda install -c omnia doxygen
conda install -c eumetsat fftw3

#~ pip install boost

### On zeolyte ###
conda install gcc

# pip install ez_setup
pip install sphinxcontrib-bibtex autodoc pyh5md

conda install gcc boost


# To instal on ZEOLITE

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
export PATH="$HOME/hydrogel/miniconda3/bin:$PATH"

conda install cmake numpy cython
#conda install cmake numpy gcc sphinx  cython hdf5  ipython
conda install -c eumetsat fftw3

conda install gcc
conda install cmake=3.6.3


OPENMPIPATH=/usr/local/programs/common/openmpi/openmpi-1.10.2/arch/x86_64-gcc_4.7.2-ofed_3.18_1
BOOSTDIR=/usr/local/programs/common/boost/boost-1.65.0/arch/x86_64-gcc_4.7.2
cmake \
-DMPI_C_LIBRARIES=$OPENMPIPATH/lib/ \
-DMPI_C_COMPILER=$OPENMPIPATH/bin/mpicc \
-DMPI_CXX_COMPILER=$OPENMPIPATH/bin/mpicxx \
-DMPI_C_INCLUDE_PATH=$OPENMPIPATH/include/ \
-DBoost_DIR=$BOOSTDIR/ \
-DBoost_INCLUDE_DIR=$BOOSTDIR/include/ ..

make -j






sudo ./cuda_9.0.176_384.81_linux.run --silent --driver --toolkit --samples --samplespath=/home/kvint/nvidia/ --run-nvidia-xconfig  --verbose

    cloog: 0.18.0-0        
    gcc:   4.8.5-7         
    gmp:   6.1.2-h6c8ec71_1
    isl:   0.12.2-0        
    mpc:   1.0.3-hec55b23_5
    mpfr:  3.1.5-h11a74b3_2



module add cmake-3.6.1 python-3.6.2-gcc


#~ To install on tarkill
git clone https://github.com/espressomd/espresso.git
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install cmake

module avail fftw
module load fftw-3.3.8
module show fftw-3.3.8
module load boost

virtualenv -p python3 venv
echo '. venv/bin/activate' > env.rc
. env.rc
pip install cython numpy ipython
cd espresso
mkdir es-build
cd es-build
cmake -DFFTW3_LIBRARIES=/software/fftw/3.3.8/lib/libfftw3.so -DFFTW3_INCLUDE_DIR=/software/fftw/3.3.8/include ..
cmake -DNUMPY_VERSION=/storage/praha1/home/kvint/miniconda3/lib/python3.6/site-packages/numpy/version.py -DNUMPY_INCLUDE_DIR=/storage/praha1/home/kvint/miniconda3/lib/python3.6/site-packages/numpy/core/include ..

#cmake ..
make -j8
