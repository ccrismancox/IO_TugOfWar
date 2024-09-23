# Testing on Ubuntu 22.04
REPDIR=`pwd`
sudo apt update
sudo apt -y upgrade


sudo apt -y install build-essential
sudo apt -y install python3-dev  git 
sudo apt -y install python-dev-is-python3
sudo apt -y install texlive-full
sudo apt -y install gfortran automake shtool libtool
sudo apt -y install python3-scipy python3-pandas python3-sympy python3-nose spyder3
sudo apt -y install subversion swig
sudo apt -y install openmpi-bin openmpi-doc libopenmpi-dev
sudo apt -y  install r-base-dev

mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4
sudo mkdir -p /etc/apt/keyrings/
wget -O ./R.asc "http://keyserver.ubuntu.com/pks/lookup?op=get&search=0xe298a3a825c0d65dfd57cbb651716619e084dab9"
gpg --no-default-keyring --keyring ./temp.gpg --import R.asc
gpg --no-default-keyring --keyring ./temp.gpg --export --output R.gpg
rm R.asc
rm temp.gpg
sudo mv R.gpg /etc/apt/keyrings
sudo echo "deb [signed-by=/etc/apt/keyrings/R.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"| sudo tee -a  /etc/apt/sources.list

sudo apt update 
sudo apt -y upgrade 

sudo apt -y install libboost-all-dev
sudo apt -y install libboost-python-dev
sudo apt -y install libboost-system-dev

sudo pip install mpi4py



cd ~
HOMEDIR=`pwd`

# This is for replication archives when this setup file needs to be separated from the git 
# ADOLC, Colpack, and IPOPT, use code stored on my bitbucket to help here
git clone https://ccrismancox@bitbucket.org/ccrismancox/pyopterf_python3.git pyopt


cd pyopt

bash setup.sh


# For first time use only
libdir=${HOMEDIR}/pyopt/Ipopt-3.12.3/lib
sudo echo $libdir$'\r' | sudo tee -a  /etc/ld.so.conf
sudo ldconfig


#use this line for project/package specific R dependencies like curl xml ssl and cmake 
sudo apt -y install  libcurl4-gnutls-dev libxml2-dev libssl-dev cmake libudunits2-dev 
cd $REPDIR
# For replications insert R package building here

echo "Ubuntu Setup Complete"

