#!/bin/bash
# Deguise Victor
# 12/05/2021

# To execute the script : ./Nettoyage_Transcriptome.bash assemblie_file R1 R2
# Don't forget R1/R2 not already create
# Refer to line 48 to know ther names. 

# DIRECTORIES

MYDIR=/mnt/c/Users/Victor/Desktop/Cours/Master/S2/Stage/Nettoyage_Transcriptomes
ASSEMBLIES_DIR=/mnt/c/Users/Victor/Desktop/Cours/Master/S2/Stage/Nettoyage_Transcriptomes/RNA-Seq_Proteogam_assemblies
FASTQ_DIR=/mnt/c/Users/Victor/Desktop/Cours/Master/S2/Stage/Nettoyage_Transcriptomes/RNA-Seq_Proteogam_Fastq
TRINITY=~/Stage/trinityrnaseq-v2.12.0/util
RESULTS=~/Stage/Results

# INITIALIZATION

FILE_GZ=$(ls $ASSEMBLIES_DIR/*.gz) 
# Check if there are files gz format
if [ ! -z "$FILE_GZ" ] ; #If variable not empty, then unzip files .gz
then
	gunzip -n $ASSEMBLIES_DIR/*.gz # Unzip all .gz files
fi


# Rename Headers for all assemblies fasta

while true; do
	read -p "Do you want to rename headers for assemblies fasta ?" yn 
	case $yn in
		[Yy]* ) for i in $(ls $ASSEMBLIES_DIR); do 
			sed -i "s/\(^>\).* \(.*_DN[0-9]*_c[0-9]*_g[0-9]*_i[0-9]*\).*/\1\2/g" $ASSEMBLIES_DIR/${i}
			done; break;;
		[Nn]* ) break;;
		* ) echo "Please answer yes or no.";;
	esac
done


# Concat Fastq
for i in $(ls $FASTQ_DIR); do
	for j in $(echo ${i} | cut -d '_' -f4 | cut -d '.' -f1 ); do # cut permet d'extraire des caractères, ici on garde le R1 ou R2.
		if [ ! -f $FASTQ_DIR/${i::4}${j}.fastq ]; # si le fichier n'existe pas alors on le créer.
			then
				echo "### Concatening Files" 
				zcat $FASTQ_DIR/${i::-12}*${j}.fastq.gz >>$FASTQ_DIR/${i::4}${j}.fastq # i::-12 pour enlever les 12 derniers caractères
				# Le fichier prendra les 4 premiers caractères du fastq lu, le R1 ou R2 puis .fastq
			fi
	done
	#rm $FASTQ_DIR/${i} # suppression du Fastq compressé
done


### CHOIX DES INPUTS ###

my_transcript=$1
R1=$2
R2=$3

# Abundance estimation 

$TRINITY/align_and_estimate_abundance.pl \
--transcripts ${my_transcript} --seqType fq \
--left ${R1} --right ${R2} --SS_lib_type RF --est_method RSEM \
--aln_method bowtie --prep_reference --output_dir rsem_outdir/${my_transcript::8} --trinity_mode \
--thread_count 20		# On fait ${my_transcript::X} pour prendre que les X 1er caractères

echo "### First step end"

# Abundance to matrix 

$TRINITY/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map none \
--out_prefix RSEM $MYDIR/rsem_outdir/${my_transcript::8}/RSEM.isoforms.results

echo "### Second step end "

# Filtering low expr highest iso 

$TRINITY/filter_low_expr_transcripts.pl \
--matrix $MYDIR/RSEM.isoform.TPM.not_cross_norm \
--transcripts ${my_transcript} \
--highest_iso_only --trinity_mode >$RESULTS/highest_iso_${my_transcript}

echo "### Finish !!"

# On fait de la place 

rm $MYDIR/RSEM.*
rm $ASSEMBLIES_DIR/${my_transcript}.*
rm -r $MYDIR/rsem_outdir/*
