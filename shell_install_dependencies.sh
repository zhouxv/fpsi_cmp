#! /bin/bash
set -e

# This script installs all necessary dependencies for the project.

install_volepsi() {
    rm -rf volepsi

    printf "################## Cloning volepsi repository   ###################\n\n"
    git clone https://github.com/ladnir/volepsi.git
    cd volepsi
    git checkout ed943f5f814591cdf864777c73b7bc9e7526c1a8

    printf "################## Building volepsi             ###################\n\n"
    sed -i '37a\ -DENABLE_FOLEAGE=ON' thirdparty/getLibOTe.cmake
    sed -i '61c\ set(libOTe_options silentot silent_vole circuits foleage)' cmake/findDependancies.cmake

    # Download boost 1.86.0, because the automatic download source of the library is too slow
    mkdir -p out && cd out
    curl -fL --retry 3 \
    -o boost_1_86_0.tar.bz2 \
    'https://sourceforge.net/projects/boost/files/boost/1.86.0/boost_1_86_0.tar.bz2/download'
    cd ..

    python3 build.py -DVOLE_PSI_ENABLE_BOOST=ON
    python3 build.py --install=../../install/volepsi
    cp ./out/build/linux/volePSI/config.h ../../install/volepsi/include/volePSI/config.h

    cd ..
}

mkdir -p thirdparty && cd thirdparty

install_volepsi
