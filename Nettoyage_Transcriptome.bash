#!/bin/bash
# Deguise Victor
# 12/05/2021


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
# Enlever le "#" devant sed si besoin de renommer les Headers. 

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
	for j in $(echo ${i} | cut -d '_' -f4 | cut -d '.' -f1 ); do
		if [ ! -f $FASTQ_DIR/${i::4}${j}.fastq ];
			then
				echo "### Concatening Files" 
				zcat $FASTQ_DIR/${i::-12}*${j}.fastq.gz >>$FASTQ_DIR/${i::4}${j}.fastq
			fi
	done
	#rm ${i} # suppression du Fastq compressé, sauf si on veut garder le fichier ?
done

# Abundance estimation 

### CHOIX DES INPUTS ###

my_transcript=$ASSEMBLIES_DIR/GPulex-f_GHCP01.1.fsa_nt
R1=$FASTQ_DIR/C4F_R1.fastq
R2=$FASTQ_DIR/C4F_R2.fastq

#mkdir -p $MYDIR/rsem_outdir/Male
#mkdir -p $MYDIR/rsem_outdir/Femelle

$TRINITY/align_and_estimate_abundance.pl \
--transcripts ${my_transcript} --seqType fq \
--left ${R1} --right ${R2} --SS_lib_type RF --est_method RSEM \
--aln_method bowtie --prep_reference --output_dir rsem_outdir/${my_transcript::8} --trinity_mode \
--thread_count 20		# On fait /${my_transcript::X} pour prendre que les X 1er caractères, ici "GPulex-f"

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
