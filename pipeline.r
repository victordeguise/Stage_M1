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


## Pour avoir acc?s a la commande makeblastdb et blast, il faut situer l'emplacement du dossier qui contient l'executable
## A t?l?charger ici : https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ (en fonction du syst?me d'exploitation)

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

list_cyp_Hya <- c("CYP18H1", "CYP302A1", "CYP306D1", "CYP307A2", "CYP314A1", "CYP315A1")
fasta_cyp_Hya <- vector(mode="list", length=6)
names(fasta_cyp_Hya) <- c("CYP18H1", "CYP302A1", "CYP306D1", "CYP307A2", "CYP314A1", "CYP315A1")

for ( i in 1:length(list_cyp_Hya) ){
  
  fasta_cyp_Hya[[i]] <- readAAStringSet(paste("Genes_Biosynthesis/",list_cyp_Hya[i],"_Hya.fasta", sep =""), format = "fasta")
  
}

### Blast contre le transcriptome femelle/male ###

cyp_blast_f <- vector(mode="list", length=6)
cyp_blast_m <- vector(mode="list", length=6)
names(cyp_blast_f) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")
names(cyp_blast_m) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")


for ( i in 1:length(list_cyp_Hya) ){
  
  cyp_blast_f[[i]] <- predict(blast_f, fasta_cyp_Hya[[i]])
  cyp_blast_m[[i]] <- predict(blast_m, fasta_cyp_Hya[[i]])
  
}

### Top 1 r?sultats Blast femelle/male ###

top_cyp_blast_f <- vector(mode="list", length=6)
top_cyp_blast_m <- vector(mode="list", length=6)
names(top_cyp_blast_f) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")
names(top_cyp_blast_m) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")

for ( i in 1:length(list_cyp_Hya) ){
  
  top_cyp_blast_f[[i]] <- cyp_blast_f[[i]][1,]
  top_cyp_blast_m[[i]] <- cyp_blast_m[[i]][1,]
  
}

### Assignation du fasta ? partir de l'id du transcriptome ###

library(seqinr)
# install.packages(seqinr)

transcriptome_f <- read.fasta("Donnees/Transcriptomes/highest_iso_GFBF_GHCZ01.1.fsa_nt")
transcriptome_m <- read.fasta("Donnees/Transcriptomes/highest_iso_GFBM_GHDA01.1.fsa_nt")

fasta_cyp_blast_f <- vector(mode="list", length=6)
fasta_cyp_blast_m <- vector(mode="list", length=6)
names(fasta_cyp_blast_f) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")
names(fasta_cyp_blast_m) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")

for ( i in 1:length(list_cyp_Hya) ){
  
  fasta_cyp_blast_f[[i]] <- transcriptome_f[top_cyp_blast_f[[i]]$SubjectID]
  fasta_cyp_blast_m[[i]] <- transcriptome_m[top_cyp_blast_m[[i]]$SubjectID]
  
}


### Recherche d'ORFs chez le male et la femelle et ?criture dans un r?pertoire ###

#install.packages("LncFinder")
library(LncFinder)

list_cyp <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")
orf_dna_f <- vector(mode="list", length=6)
orf_dna_m <- vector(mode="list", length=6)
orf_Forward_aa_cyp_f = list()
orf_Forward_aa_cyp_m = list()
names(orf_dna_f) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")
names(orf_dna_m) <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")

for ( i in 1:length(list_cyp) ){
  
  orf_dna_f[[i]] <- find_orfs(fasta_cyp_blast_f[[i]],reverse.strand = FALSE)
  orf_dna_m[[i]] <- find_orfs(fasta_cyp_blast_m[[i]],reverse.strand = FALSE)
  orf_Forward_aa_cyp_f[[i]] <- paste(translate(s2c(orf_dna_f[[i]]$ORF.Max.Seq)), collapse = "")
  orf_Forward_aa_cyp_m[[i]] <- paste(translate(s2c(orf_dna_m[[i]]$ORF.Max.Seq)), collapse = "")
  # Ecriture des r?sultats dans le dossier R_Orf/
  write.fasta(orf_Forward_aa_cyp_f[[i]], names = top_cyp_blast_f[[i]]$SubjectID, 
              file.out = paste("R_Orf/orf_Forward_",list_cyp[i],"_f.fasta", sep=""), open = "w")
  write.fasta(orf_Forward_aa_cyp_m[[i]], names = top_cyp_blast_m[[i]]$SubjectID, 
              file.out = paste("R_Orf/orf_Forward_",list_cyp[i],"_m.fasta", sep=""), open = "w")
}


### Alignement multiple pour la phylog?nie ###

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
library(msa)
library(ape)

## Lire dux fichier fasta pour mettre notre ORF a la fin du fichier d'alignement

alignement <- readAAStringSet(c("Phylogeny/Arthropod_phylogeny.fasta","R_Orf/orf_forward_cyp18_f.fasta" ,
                                "R_Orf/orf_forward_cyp302_f.fasta","R_Orf/orf_forward_cyp306_f.fasta", 
                                "R_Orf/orf_forward_cyp306_m.fasta", "R_Orf/orf_forward_cyp307_f.fasta", 
                                "R_Orf/orf_forward_cyp314_f.fasta", "R_Orf/orf_forward_cyp315_f.fasta"))

multiple_alignement <- msaClustalW(alignement, type = "protein")

## Si on veut garder l'alignement multiple : 
#write.phylip(multiple_alignement, "Phylogeny/Multiple_alignement.txt")

#msaPrettyPrint(multiple_alignement, output = "pdf")



### Phylogénie ###

multiple_alignement <- msaConvert(multiple_alignement)
d <- dist.alignment(multiple_alignement, "identity")
tree <- nj(d)
tree <- makeLabel(tree, space = "")
plot.phylo(tree,type = "phylogram", main="Phylogenetic Tree", use.edge.length = FALSE, font = 2 )




### R?cup?rer les 10 top blast Hits du cyp306 puis rechercher les ORFs  

top_cyp306_blast_f <- cyp306_blast_f[1:10,]
top_cyp306_blast_m <- cyp306_blast_m[1:10,]
fasta_cyp306_blast_f <- transcriptome_f[top_cyp306_blast_f$SubjectID]
fasta_cyp306_blast_m <- transcriptome_m[top_cyp306_blast_m$SubjectID]




for ( i in 1:10 ) { 

  orf_dna_cyp306_f[[i]] <- find_orfs(fasta_cyp306_blast_f[[i]], reverse.strand = TRUE)
  orf_Forward_aa_cyp306_f[[i]] <-  paste(translate(s2c(orf_dna_cyp306_f[[i]]$ORF.Forward$ORF.Max.Seq)), collapse = "")
  write.fasta(orf_Forward_aa_cyp306_f[[i]], names = top_cyp306_blast_f[i,]$SubjectID, 
              file.out = paste("R_Orf/orf_forward_cyp306_f", i,".fasta", sep=""), open = "w")
  orf_Reverse_aa_cyp306_f[[i]] <-  paste(translate(s2c(orf_dna_cyp306_f[[i]]$ORF.Reverse$ORF.Max.Seq)), collapse = "")
  write.fasta(orf_Reverse_aa_cyp306_f[[i]], names = top_cyp306_blast_f[i,]$SubjectID, 
              file.out = paste("R_Orf/orf_Reverse_cyp306_f", i,".fasta", sep=""), open = "w")
}


for ( i in 1:10 ) { 
  
  orf_dna_cyp306_m[[i]] <- find_orfs(fasta_cyp306_blast_m[[i]], reverse.strand = TRUE)
  orf_Forward_aa_cyp306_m[[i]] <-  paste(translate(s2c(orf_dna_cyp306_m[[i]]$ORF.Forward$ORF.Max.Seq)), collapse = "")
  write.fasta(orf_Forward_aa_cyp306_m[[i]], names = top_cyp306_blast_m[i,]$SubjectID, 
              file.out = paste("R_Orf/orf_forward_cyp306_m", i,".fasta", sep=""), open = "w")
  orf_Reverse_aa_cyp306_m[[i]] <-  paste(translate(s2c(orf_dna_cyp306_m[[i]]$ORF.Reverse$ORF.Max.Seq)), collapse = "")
  write.fasta(orf_Reverse_aa_cyp306_m[[i]], names = top_cyp306_blast_m[i,]$SubjectID, 
              file.out = paste("R_Orf/orf_Reverse_cyp306_m", i,".fasta", sep=""), open = "w")
}








