ssh -i /drives/d/ssh/denbi.key centos@193.196.20.111 # mobaxterm

chandler-1

installed programs

R
miniconda
tree
byobu
bzip2
mothur

# create and attach a hard drive under horizon dashboard, 
# then find your unformatted volume (check size), here "vdb"
lsblk
# and format volume with format ext4
mkfs.ext4 /dev/vdb
# make it accessible for your vm by creating a mountpoint
mkdir -p /mnt/data
# grant yourself admin rights 
chmod 777 /mnt/data/  # 777 explained: https://en.wikipedia.org/wiki/Chmod#Numeric_example
# link the formatted volumne "vdb" to the folder "mnt/data"
mount /dev/vdb /mnt/data
# check if it is there using
df -h
# unmount with
umount /dev/vdb

# install programs on centOS
sudo yum install bzip
sudo wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
mkdir /data/downloads
mv Miniconda3-latest-Linux-x86_64 /data/downloads
sh /data/downloads/Miniconda3-latest-Linux-x86_64.sh 
# make sure to restart terminal afterwards

# first, make sure your conda channels are configured correctly
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install R and dada2 into conda environment

conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8
cd /data
mkdir projects
mkdir programs

# install mothur 1.41.0
cd /data/downloads
wget https://github.com/mothur/mothur/releases/download/v1.41.0.pre-release/Mothur.linux_64.zip
unzip Mothur.linux_64.zip
rm -r __MACOSX
mv ./downloads/mothur/ ./programs/
# run mothur
/data/programs/mothur/mothur

# setup projects structure for raw reads
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/all_runs
mkdir -p /data/projects/glyphosate/reads/raw_reads_metagenome



# copy data from bio-48 to chandler-1, being on chandler-1

# from local to chandler-1
scp -i /drives/d/ssh/denbi.key /drives/d/backup_sequenzen/glyphosate/amplicon*/*.gz centos@193.196.20.111:/data/projects/glyphosate/reads/raw_reads_16S/all_runs


cp -r /data/projects
scp -r janssen@10.11.20.48:/data/projects .
scp -i /drives/d/ssh/denbi.key -r janssen@10.11.20.48:/data/projects/* /data/projects/ 

# running on phy-2 (10.11.80.2?)
Detach the screen using CTRL+A d