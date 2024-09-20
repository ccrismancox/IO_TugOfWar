# Testing on Ubuntu 20.04.4 Focal 
REPDIR=`pwd`
cd ~
HOMEDIR=`pwd`
sudo apt update
sudo apt -y upgrade


sudo apt -y install build-essential
sudo apt -y install python3-dev  git 
sudo apt -y install python3dev-is-python
sudo apt -y install texlive-full
sudo apt -y install gfortran automake shtool libtool
sudo apt -y install python3-matplotlib python3-scipy python3-pandas python3-sympy python3-nose spyder3
sudo apt -y install subversion swig
sudo apt -y install openmpi-bin openmpi-doc libopenmpi-dev
sudo apt -y  install r-base-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4


sudo apt -y install libboost-all-dev
sudo apt -y install libboost-python-dev
sudo apt -y install libboost-system-dev

sudo pip install rpy2   
#==2.8.6 #no idea if this is still required; under 18.04 and python 2.7 not using this version broke things
sudo pip install mpi4py

# This is for replication archives when this setup file needs to be separated from the git 
# ADOLC, Colpack, and IPOPT, use code stored on my bitbucket to help here
git clone https://ccrismancox@bitbucket.org/ccrismancox/pyopterf_python3.git pyopt


cd pyopt

bash setup.sh


# For first time use only
libdir=${HOMEDIR}/pyopt/Ipopt-3.12.3/lib
sudo echo $libdir$'\r' | sudo tee -a  /etc/ld.so.conf
sudo ldconfig

sudo apt -y install libgfortran3 #where does this go?
#use this line for project/package specific R dependencies like curl xml ssl and cmake 
sudo apt -y install  libcurl4-gnutls-dev libxml2-dev libssl-dev cmake 
cd $REPDIR
# For replications insert R package building here

echo "Ubuntu Setup Complete"

