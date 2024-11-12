#!/bin/bash

echo "choose the group you want to process and press [ENTER] :"
ls /ChIP_seq/
read GROUP
cd /ChIP_seq/$GROUP 2> /dev/null

##############################################################################
until [ "$PWD" == "/ChIP_seq/$GROUP" ]; do
#the space in the bracket (at the begining and at the end) is very important##

echo "Please choose a real group"
ls /ChIP_seq/
read GROUP
cd /ChIP_seq/$GROUP 2> /dev/null;
done

echo "choose the sex you want to process and press [ENTER] :"
ls /ChIP_seq/$GROUP
read HISTONE
cd /ChIP_seq/$GROUP/$HISTONE 2> /dev/null
until [ "$PWD" == "/ChIP_seq/$GROUP/$HISTONE" ]; do
echo "Please choose a real sex mark"
ls /ChIP_seq/
read HISTONE
cd /ChIP_seq/$GROUP/$HISTONE 2> /dev/null;
done

echo "choose the subgroup you want to process and press [ENTER] :"
ls /ChIP_seq/$GROUP/$HISTONE
read SUBGROUP
cd /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP 2> /dev/null
until [ "$PWD" == "/ChIP_seq/$GROUP/$HISTONE/$SUBGROUP" ]; do
echo "Please choose a real subgroup"
ls /ChIP_seq/$GROUP/$HISTONE
read SUBGROUP
cd /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP 2> /dev/null;
done

cd /ChIP_seq/
sed "s/GROUPGOESHERE/"$GROUP"/" ./ChIP-seq_script.sh > /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_1.sh

sed "s/SUBgroupGOESHERE/"$SUBGROUP"/" /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_1.sh > /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_2.sh

sed "s/HISTONEGOESHERE/"$HISTONE"/" /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr_2.sh > /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP/precorr.sh

cd /ChIP_seq/$GROUP/$HISTONE/$SUBGROUP

rm ./precorr_1.sh ./precorr_2.sh
chmod u+x ./precorr.sh

sbatch -p normal_q --account=chipseq -t 0-12:00:00 --mem=20GB ./precorr.sh
