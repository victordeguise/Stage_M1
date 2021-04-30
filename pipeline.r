#Victor Deguise
#27/04/2021


### INITIALISATION ###

setwd(dir= "C:/Users/Victor/Desktop/Cours/Master/S2/Stage/Projet/")
dir= "C:/Users/Victor/Desktop/Cours/Master/S2/Stage/Projet/"


# To install blast :
#Install the package "devtools"
#install.packages("devtools")
#Install the Bioconductor package Biostrings either using devtools::install_bioc("Biostrings")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#Install rBlast from GitHub using devtools::install_github("mhahsler/rBLAST").
# Pour afficher l'alignement multiple dans un pdf il faut ce package
#install.packages('tinytex')




Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:\\Users\\Victor\\Desktop\\Cours\\Master\\S2\\Stage\\Projet\\ncbi-blast-2.11.0+\\bin",
                        sep= .Platform$path.sep))
library(rBLAST)

### Cr?ation des bases de donn?es avec nos transcriptomes femelle et male ###

makeblastdb(file="Donnees/Transcriptomes/highest_iso_GFBF_GHCZ01.1.fsa_nt", dbtype = "nucl", args = c("-out GfossB_f -title GfossB_f -parse_seqids"))
makeblastdb(file="Donnees/Transcriptomes/highest_iso_GFBM_GHDA01.1.fsa_nt", dbtype = "nucl", args=c("-out GfossB_m -title GfossB_m -parse_seqids"))

### load a blast database ###

blast_f <- blast(db="GfossB_f", type= "tblastn")
blast_m <- blast(db="GfossB_m", type= "tblastn")

### Affiliation des fichiers fasta ###

cyp302 <- readAAStringSet("Genes_Biosynthesis/CYP302A1_Hya.fasta", format = "fasta")
cyp306 <- readAAStringSet("Genes_Biosynthesis/CYP306D1_Hya.fasta", format = "fasta")
cyp307 <- readAAStringSet("Genes_Biosynthesis/CYP307A2_Hya.fasta", format = "fasta")
cyp314 <- readAAStringSet("Genes_Biosynthesis/CYP314A1_Hya.fasta", format = "fasta")
cyp315 <- readAAStringSet("Genes_Biosynthesis/CYP315A1_Hya.fasta", format = "fasta")

### Blast contre le transcriptome femelle ###

cyp302_blast_f <- predict(blast_f, cyp302)
cyp306_blast_f <- predict(blast_f, cyp306)
cyp307_blast_f <- predict(blast_f, cyp307)
cyp314_blast_f <- predict(blast_f, cyp314)
cyp315_blast_f <- predict(blast_f, cyp315)


### Blast contre le transcriptome male ###

cyp302_blast_m <- predict(blast_m, cyp302)
cyp306_blast_m <- predict(blast_m, cyp306)
cyp307_blast_m <- predict(blast_m, cyp307)
cyp314_blast_m <- predict(blast_m, cyp314)
cyp315_blast_m <- predict(blast_m, cyp315)


### Top 1 r?sultats Blast femelle ###

top_cyp302_blast_f <- cyp302_blast_f[1,]
top_cyp306_blast_f <- cyp306_blast_f[1,]
top_cyp307_blast_f <- cyp307_blast_f[1,]
top_cyp314_blast_f <- cyp314_blast_f[1,]
top_cyp315_blast_f <- cyp315_blast_f[1,]


### Top 1 r?sultats Blast male ###

top_cyp302_blast_m <- cyp302_blast_m[1,]
top_cyp306_blast_m <- cyp306_blast_m[1,]
top_cyp307_blast_m <- cyp307_blast_m[1,]
top_cyp314_blast_m <- cyp314_blast_m[1,]
top_cyp315_blast_m <- cyp315_blast_m[1,]


### Assignation du fasta ? partir de l'id du transcriptome ###

### pour la femelle ###

library(seqinr)
# install.packages(seqinr)
transcriptome_f <- read.fasta("Donnees/Transcriptomes/highest_iso_GFBF_GHCZ01.1.fsa_nt")

fasta_cyp302_blast_f <- transcriptome_f[top_cyp302_blast_f$SubjectID]
fasta_cyp306_blast_f <- transcriptome_f[top_cyp306_blast_f$SubjectID]
fasta_cyp307_blast_f <- transcriptome_f[top_cyp307_blast_f$SubjectID]
fasta_cyp314_blast_f <- transcriptome_f[top_cyp314_blast_f$SubjectID]
fasta_cyp315_blast_f <- transcriptome_f[top_cyp315_blast_f$SubjectID]


### pour le male ### 
transcriptome_m <- read.fasta("Donnees/Transcriptomes/highest_iso_GFBM_GHDA01.1.fsa_nt")

fasta_cyp302_blast_m <- transcriptome_m[top_cyp302_blast_m$SubjectID]
fasta_cyp306_blast_m <- transcriptome_m[top_cyp306_blast_m$SubjectID]
fasta_cyp307_blast_m <- transcriptome_m[top_cyp307_blast_m$SubjectID]
fasta_cyp314_blast_m <- transcriptome_m[top_cyp314_blast_m$SubjectID]
fasta_cyp315_blast_m <- transcriptome_m[top_cyp315_blast_m$SubjectID]


### Recherche d'ORFs ###

#install.packages("LncFinder")
library(LncFinder)

# ORF chez la femelle

orf_dna_cyp302_f <-find_orfs(fasta_cyp302_blast_f,reverse.strand = TRUE)
orf_Forward_aa_cyp302_f <- paste(translate(s2c(orf_dna_cyp302_f$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp302_f <- paste(translate(s2c(orf_dna_cyp302_f$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp306_f <-find_orfs(fasta_cyp306_blast_f,reverse.strand = TRUE )
orf_Forward_aa_cyp306_f <- paste(translate(s2c(orf_dna_cyp306_f$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp306_f <- paste(translate(s2c(orf_dna_cyp306_f$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp307_f <-find_orfs(fasta_cyp307_blast_f,reverse.strand = TRUE )
orf_Forward_aa_cyp307_f <- paste(translate(s2c(orf_dna_cyp307_f$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp307_f <- paste(translate(s2c(orf_dna_cyp307_f$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp314_f <-find_orfs(fasta_cyp314_blast_f,reverse.strand = TRUE )
orf_Forward_aa_cyp314_f <- paste(translate(s2c(orf_dna_cyp314_f$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp314_f <- paste(translate(s2c(orf_dna_cyp314_f$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp315_f <-find_orfs(fasta_cyp315_blast_f,reverse.strand = TRUE )
orf_Forward_aa_cyp315_f <- paste(translate(s2c(orf_dna_cyp315_f$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp315_f <- paste(translate(s2c(orf_dna_cyp315_f$ORF.Reverse$ORF.Max.Seq)), collapse = "")

# ORF chez le male

orf_dna_cyp302_m <-find_orfs(fasta_cyp302_blast_m,reverse.strand = TRUE )
orf_Forward_aa_cyp302_m <- paste(translate(s2c(orf_dna_cyp302_m$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp302_m <- paste(translate(s2c(orf_dna_cyp302_m$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp306_m <-find_orfs(fasta_cyp306_blast_m,reverse.strand = TRUE )
orf_Forward_aa_cyp306_m <- paste(translate(s2c(orf_dna_cyp306_m$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp306_m <- paste(translate(s2c(orf_dna_cyp306_m$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp307_m <-find_orfs(fasta_cyp307_blast_m,reverse.strand = TRUE )
orf_Forward_aa_cyp307_m <- paste(translate(s2c(orf_dna_cyp307_m$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp307_m <- paste(translate(s2c(orf_dna_cyp307_m$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp314_m <-find_orfs(fasta_cyp314_blast_m,reverse.strand = TRUE )
orf_Forward_aa_cyp314_m <- paste(translate(s2c(orf_dna_cyp314_m$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp314_m <- paste(translate(s2c(orf_dna_cyp314_m$ORF.Reverse$ORF.Max.Seq)), collapse = "")

orf_dna_cyp315_m <-find_orfs(fasta_cyp315_blast_m,reverse.strand = TRUE )
orf_Forward_aa_cyp315_m <- paste(translate(s2c(orf_dna_cyp315_m$ORF.Forward$ORF.Max.Seq)), collapse = "")
orf_Reverse_aa_cyp315_m <- paste(translate(s2c(orf_dna_cyp315_m$ORF.Reverse$ORF.Max.Seq)), collapse = "")

### Ecriture des r?sultats dans le dossier R_Orf ###

### pour la femelle
### Forward strand

write.fasta(orf_Forward_aa_cyp302_f, names = top_cyp302_blast_f$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp302_f.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp306_f, names = top_cyp306_blast_f$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp306_f.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp307_f, names = top_cyp307_blast_f$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp307_f.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp314_f, names = top_cyp314_blast_f$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp314_f.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp315_f, names = top_cyp315_blast_f$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp315_f.fasta", open = "w")

### Reverse strand

write.fasta(orf_Reverse_aa_cyp302_f, names = top_cyp302_blast_f$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp302_f.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp306_f, names = top_cyp306_blast_f$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp306_f.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp307_f, names = top_cyp307_blast_f$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp307_f.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp314_f, names = top_cyp314_blast_f$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp314_f.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp315_f, names = top_cyp315_blast_f$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp315_f.fasta", open = "w")



### pour le male
## Forward strand 

write.fasta(orf_Forward_aa_cyp302_m, names = top_cyp302_blast_m$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp302_m.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp306_m, names = top_cyp306_blast_m$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp306_m.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp307_m, names = top_cyp307_blast_m$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp307_m.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp314_m, names = top_cyp314_blast_m$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp314_m.fasta", open = "w")
write.fasta(orf_Forward_aa_cyp315_m, names = top_cyp315_blast_m$SubjectID, 
            file.out = "R_Orf/orf_forward_cyp315_m.fasta", open = "w")

### Reverse strand

write.fasta(orf_Reverse_aa_cyp302_m, names = top_cyp302_blast_m$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp302_m.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp306_m, names = top_cyp306_blast_m$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp306_m.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp307_m, names = top_cyp307_blast_m$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp307_m.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp314_m, names = top_cyp314_blast_m$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp314_m.fasta", open = "w")
write.fasta(orf_Reverse_aa_cyp315_m, names = top_cyp315_blast_m$SubjectID, 
            file.out = "R_Orf/orf_Reverse_cyp315_m.fasta", open = "w")


### Alignement multiple pour la phylog?nie ###

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
library(msa)
library(ape)

multiple_alignement <- msaMuscle(readAAStringSet("Phylogeny/First_Phylogeny_Female.fasta"))
write.phylip(multiple_alignement, "Phylogeny/Multiple_alignement.txt")


msaPrettyPrint(multiple_alignement, output = "pdf")

multiple_alignement <- msaConvert(multiple_alignement)

d <- dist.alignment(multiple_alignement, "identity")
plot(nj(d), main="Phylogenetic Tree") 



### Phylogénie ###

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggtree)







