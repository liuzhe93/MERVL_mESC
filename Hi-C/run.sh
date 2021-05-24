1.install Juicer
mkdir ./opt cd opt
git clone https://github.com/theaidenlab/juicer.git
ln -s juicer/CPU scripts cd scripts/common
wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
# index reference (Here, human)
mkdir references cd references cp <path>/Homo_sapiens_assembly38.fasta  
# OR ln -s <path>/Homo_sapiens_assembly38.fasta bwa 
index  Homo_sapiens_assembly38.fasta
index  mm9.fasta
cd .. 
# Ìí¼Órestriction enzyme site (Here, MboI)
mkdir restriction_sites
cd restriction_sites/ python2 ../juicer/misc/generate_site_positions.py MboI  mm9_MboI ../references/mm9.fasta 
awk 'BEGIN{OFS="\t"}{print $1, $NF}' mm9_MboI.txt > mm9.chrom.sizes
cd ..
2. Import Hi-C raw data 
mkdir fastq
cd fastq
# put your data in this direction
# note that: the format of data should be:  _R1.fastq,  _R2.fastq
cd .. 
for i in ¡°treatment1 treatment2 treatment3¡±
3.Hi-C data analysis using Juicer software
do 
    <path>/juicer/scripts/juicer.sh -y <path>/juicer/restriction_sites/file/hg19_MboI.txt -p /home/liuzhe/Data/genome/chrom.sizes/hg19.chrom.sizes -d <path>juicer/fastq/$i/-z /home/liuzhe/YanLab/Project3_Hi-C/SmallSampleTest_juicer/juicer/references/mm9.fasta -D /home/liuzhe/YanLab/Project3_Hi-C/SmallSampleTest_juicer/juicer -x
done
#-y: restriction site file;-p: chrom.sizes path; -d: topDir;-D: Juicer scripts directory;-x: exclude fragment-delimited maps from hic file creation
# input: .fastq file; output: .hic file
# inter_30.txt can tell you all your commands and results information statistics. _30: filtered by using MAQ>30. 
4. [optional] Identify TADs
java -jar <path>/juicer/scripts/juicer_tools.jar arrowhead -m 2000 -r 10000 aligned/inter_30.hic aligned/res10kb/treatment.bedpe
# -m: matrix size; -r: resolution
5. [optional] Convert .hic to matrix file
Juicer_matrix.py