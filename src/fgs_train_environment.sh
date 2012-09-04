#!/bin/bash
# Prepare environment for FGS retraining 

sudo apt-get update
sudo apt-get install git
sudo apt-get install python-biopython -y
sudo chown ubuntu /mnt

genomedir=/mnt/genomes

mkdir build
cd build
git clone https://github.com/wtangiit/fgs_train.git
git clone https://github.com/wltrimbl/FGS.git

echo export PATH=\$PATH:$HOME/build/FGS:$HOME/build/fgs_train/src >> ~/.profile

mkdir $genomedir
curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz > $genomedir/all.fna.tar.gz
curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.ptt.tar.gz > $genomedir/all.ptt.tar.gz
tar xvf all.fna.tar.gz 
tar xvf all.ptt.tar.gz

ln -s genomes/Acinetobacter_baumannii_AB307_0294_uid59271/* /mnt
cd /mnt

gen_train_input.py -i NC_011595.fna -p NC_011595.ptt
# generates 
#    train_fwd_NC_011595.fna.csv
#    train_noncoding_NC_011595.fna.csv

validate_train_input.py -i train_fwd_NC_011595.fna.csv | tee codonstats.txt

fgs_train.py  -i train_fwd_NC_011595.fna.csv -n train_noncoding_NC_011595.fna.csv
# generates gene, rgene, start, start1, stop, stop1, gene.ct, rgene.ct

