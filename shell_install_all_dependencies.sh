#! /bin/bash
set -e

# This script installs all necessary dependencies for the project.

install_securejoin(){
    rm -rf secure-join

    printf "################## Cloning secure-join repository   ###################\n\n"
    git clone https://github.com/ladnir/secure-join.git
    cd secure-join
    git checkout 377ca63b9d8f4f6aede0d3a2e3d9078973a3ee10
    
    printf "################## Building secure-join             ###################\n\n"
    python3 build.py -DFETCH_SODIUM=ON -DOC_THIRDPARTY_HINT="$(pwd)/out/install/linux"

    python3 build.py --install=../../install/secure-join

    cd ..
}



install_volepsi(){
    rm -rf volepsi

    printf "################## Cloning volepsi repository   ###################\n\n"
    git clone https://github.com/ladnir/volepsi.git
    cd volepsi
    git checkout ed943f5f814591cdf864777c73b7bc9e7526c1a8

    printf "################## Building volepsi             ###################\n\n"
    python3 build.py -DVOLE_PSI_ENABLE_BOOST=ON --install=../../install/volepsi
    cp ./build/volePSI/config.h ../../install/volepsi/include/volePSI/config.h

    cd ..
}

mkdir -p thirdparty && cd thirdparty


install_securejoin
install_volepsi


if [ -f /.dockerenv ]; then
    echo "Running inside a Docker container, cleaning up thirdparty source directories"
    rm -rf thirdparty
else
    echo "Not running inside a Docker container"
fi

