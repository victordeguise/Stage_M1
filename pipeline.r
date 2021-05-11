#####################################################################
## PROJECT : G. fossarum                                           ##
## STUDIES : identification cytochrome P450                        ##
## AUTHOR : Victor Deguise                                         ##
## DATE : May 2021                                                 ##
## SCRIPT : pipeline Phylogeny                                     ##
#####################################################################
#-------------------------------------------------------------------
#  INTRODUCTORY NOTE                            
#-------------------------------------------------------------------
# This script allows to identify cytochrome P450 by phylogenitical
# analysis of the G. fossarum 

# The sessionInfo() is in the appendice part

#-------------------------------------------------------------------
#  INDEX.SUMMARY OF THE PIPELINE                   
#-------------------------------------------------------------------

#      PACKAGES & SETTINGS 
#      DIRECTORIES  

# 01.  CRÉATION DES BASES DE DONNÉES AVEC NOS TRANSCRIPTOMES
# 02.  LOAD A BLAST DATABASE
# 03.  AFFILIATION DES FICHIERS FASTA
# 04.  BLAST CONTRE LE TRANSCRIPTOME FEMELLE/MALE
# 05.  SÉLECTION TOP 1 RÉSULTATS BLAST FEMELLE/MALE
# 06.  ASSIGNATION DU FASTA À PARTIR DE L'ID DU TRANSCRIPTOME
# 07.  RECHERCHE D'ORFS ET ÉCRITURE DES FASTA
# 08.  ALIGNEMENT MULTIPLE ET PHYLOGÉNIE

#      APPENDICES

#-------------------------------------------------------------------
#  PACKAGES & SETTINGS                          
#-------------------------------------------------------------------
## Required packages
# Installs missing libraries !
list.of.packages <- c("BiocManager", "devtools", "seqinr", "LncFinder", "ape", "tools", "ggtree") #list of packages required
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #list of packages non installed
if(length(new.packages)) install.packages(new.packages, repos='https://cran.rstudio.com/') #install packages if the new.packages list is not empty
#devtools::install_bioc("Biostrings")
#devtools::install_github("mhahsler/rBLAST")
#BiocManager::install("msa")

#tools
library(rBLAST)
library(seqinr)
library(LncFinder)
library(msa)
library(ape)
library(tools)
library(ggtree)


#set environment
## Pour avoir accès a la commande makeblastdb et blast, il faut situer l'emplacement du dossier qui contient l'executable
## A télécharger ici : https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ (en fonction du système d'exploitation)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:\\Users\\Victor\\Desktop\\Cours\\Master\\S2\\Stage\\Projet\\ncbi-blast-2.11.0+\\bin",sep= .Platform$path.sep))

## Pour visualiser l'alignement multiple il faut d'abord installer miktex
# https://miktex.org/download
## Puis il faut ajouter la variable d'environnement
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Users/Victor/AppData/Local/Programs/MiKTeX/miktex/bin/x64",sep= .Platform$path.sep))

#-------------------------------------------------------------------
#  DIRECTORIES                          
#-------------------------------------------------------------------

## Working directory
wdir <- getwd()
wdir #current directory

## Input directories
Transcriptome_dir <- file.path(wdir, "Donnees/Transcriptomes/")
Genes_dir <- file.path(wdir, "Genes_Biosynthesis/")

## Output directory
ORFdir <- file.path(wdir, "R_Orf/")
PhylogenyDir <- file.path(wdir, "Phylogeny/")
Multiple_align_Dir <- file.path(wdir, "Multiple_Alignement/")



#-------------------------------------------------------------------
# 01. CRÉATION DES BASES DE DONNÉES AVEC NOS TRANSCRIPTOMES 
#-------------------------------------------------------------------

transcriptome_f <- read.fasta(file.path(Transcriptome_dir, "highest_iso_GFBF_GHCZ01.1.fsa_nt"))
transcriptome_m <- read.fasta(file.path(Transcriptome_dir, "highest_iso_GFBM_GHDA01.1.fsa_nt"))

makeblastdb(file.path(Transcriptome_dir, "highest_iso_GFBF_GHCZ01.1.fsa_nt"), dbtype = "nucl", args = c("-out GfossB_f -title GfossB_f -parse_seqids"))
makeblastdb(file.path(Transcriptome_dir, "highest_iso_GFBM_GHDA01.1.fsa_nt"), dbtype = "nucl", args=c("-out GfossB_m -title GfossB_m -parse_seqids"))

#-------------------------------------------------------------------
# 02. LOAD A BLAST DATABASE 
#-------------------------------------------------------------------

blast_f <- blast(db="GfossB_f", type= "tblastn")
blast_m <- blast(db="GfossB_m", type= "tblastn")

#-------------------------------------------------------------------
# 03. AFFILIATION DES FICHIERS FASTA 
#-------------------------------------------------------------------

list_cyp_Hya <- c("CYP18H1", "CYP302A1", "CYP306D1", "CYP307A2", "CYP314A1", "CYP315A1")
fasta_cyp_Hya <- vector(mode="list", length=6)
names(fasta_cyp_Hya) <- list_cyp_Hya

for ( i in 1:length(list_cyp_Hya) ){
  
  fasta_cyp_Hya[[i]] <- readAAStringSet(paste(Genes_dir,"/",list_cyp_Hya[i],"_Hya.fasta", sep =""), format = "fasta")
  
}

#-------------------------------------------------------------------
# 04. BLAST CONTRE LE TRANSCRIPTOME FEMELLE/MALE
#-------------------------------------------------------------------

list_cyp <- c("CYP18", "CYP302", "CYP306", "CYP307", "CYP314", "CYP315")
cyp_blast_f <- vector(mode="list", length=6)
cyp_blast_m <- vector(mode="list", length=6)
names(cyp_blast_f) <- list_cyp
names(cyp_blast_m) <- list_cyp


for ( i in 1:length(list_cyp_Hya) ){
  
  cyp_blast_f[[i]] <- predict(blast_f, fasta_cyp_Hya[[i]])
  cyp_blast_m[[i]] <- predict(blast_m, fasta_cyp_Hya[[i]])
  
}

#-------------------------------------------------------------------
# 05. SÉLECTION TOP 1 RÉSULTATS BLAST FEMELLE/MALE
#-------------------------------------------------------------------

top_cyp_blast_f <- vector(mode="list", length=6)
top_cyp_blast_m <- vector(mode="list", length=6)
names(top_cyp_blast_f) <- list_cyp
names(top_cyp_blast_m) <- list_cyp

for ( i in 1:length(list_cyp_Hya) ){
  
  top_cyp_blast_f[[i]] <- cyp_blast_f[[i]][1,]
  top_cyp_blast_m[[i]] <- cyp_blast_m[[i]][1,]
  
}

#-------------------------------------------------------------------
# 06. ASSIGNATION DU FASTA À PARTIR DE L'ID DU TRANSCRIPTOME
#-------------------------------------------------------------------

fasta_cyp_blast_f <- vector(mode="list", length=6)
fasta_cyp_blast_m <- vector(mode="list", length=6)
names(fasta_cyp_blast_f) <- list_cyp
names(fasta_cyp_blast_m) <- list_cyp

for ( i in 1:length(list_cyp_Hya) ){
  
  fasta_cyp_blast_f[[i]] <- transcriptome_f[top_cyp_blast_f[[i]]$SubjectID]
  fasta_cyp_blast_m[[i]] <- transcriptome_m[top_cyp_blast_m[[i]]$SubjectID]
  
}


#-------------------------------------------------------------------
# 07. RECHERCHE D'ORFS ET ÉCRITURE DES FASTA
#-------------------------------------------------------------------

orf_dna_f <- vector(mode="list", length=6)
orf_dna_m <- vector(mode="list", length=6)
orf_Forward_aa_cyp_f = vector(mode="list", length=6)
orf_Forward_aa_cyp_m = vector(mode="list", length=6)
names(orf_dna_f) <- list_cyp
names(orf_dna_m) <- list_cyp
names(orf_Forward_aa_cyp_f) <- list_cyp
names(orf_Forward_aa_cyp_m) <- list_cyp


for ( i in 1:length(list_cyp) ){
  
  orf_dna_f[[i]] <- find_orfs(fasta_cyp_blast_f[[i]],reverse.strand = FALSE)
  orf_dna_m[[i]] <- find_orfs(fasta_cyp_blast_m[[i]],reverse.strand = FALSE)
  orf_Forward_aa_cyp_f[[i]] <- paste(translate(s2c(orf_dna_f[[i]]$ORF.Max.Seq)), collapse = "")
  orf_Forward_aa_cyp_m[[i]] <- paste(translate(s2c(orf_dna_m[[i]]$ORF.Max.Seq)), collapse = "")
  # Ecriture des résultats dans le dossier R_Orf/
  write.fasta(orf_Forward_aa_cyp_f[[i]], names = top_cyp_blast_f[[i]]$SubjectID, 
              file.out = file.path(ORFdir, paste("orf_",list_cyp[i],"_f.fasta", sep="")), open = "w")
  write.fasta(orf_Forward_aa_cyp_m[[i]], names = top_cyp_blast_m[[i]]$SubjectID, 
              file.out = file.path(ORFdir, paste("orf_",list_cyp[i],"_m.fasta", sep="")), open = "w")
}


#-------------------------------------------------------------------
# 08. ALIGNEMENT MULTIPLE ET PHYLOGÉNIE 
#-------------------------------------------------------------------

list_orf <- list.files("R_Orf/")
list_cyp_phylo <- c("CYP18_F_", "CYP18_M_","CYP302_F_","CYP302_M_", "CYP306_F_", "CYP306_M_","CYP307_F_", 
                    "CYP307_M_","CYP314_F_", "CYP314_M_","CYP315_F_","CYP315_M_")

top_cyp_blast = rbind(top_cyp_blast_f$CYP18$SubjectID,top_cyp_blast_m$CYP18$SubjectID,top_cyp_blast_f$CYP302$SubjectID,
                      top_cyp_blast_m$CYP302$SubjectID,top_cyp_blast_f$CYP306$SubjectID,top_cyp_blast_m$CYP306$SubjectID,
                      top_cyp_blast_f$CYP307$SubjectID,top_cyp_blast_m$CYP307$SubjectID,top_cyp_blast_f$CYP314$SubjectID,
                      top_cyp_blast_m$CYP314$SubjectID,top_cyp_blast_f$CYP315$SubjectID,top_cyp_blast_m$CYP315$SubjectID)


for ( i in 1:length(list_orf)) {
  sequence <- readAAStringSet(c("Phylogeny/Phylogeny_Fasta_Arthropods.fasta",file.path(ORFdir, list_orf[i])) )
  multiple_alignement <- msa(sequence, method = "Muscle", type = "protein")
  #msaPrettyPrint(multiple_alignement, file = file.path(Multiple_align_Dir, paste(list_orf[i], ".tex", sep = "")) , output = "tex")
  #texi2pdf(file.path(Multiple_align_Dir, paste(list_orf[i], ".tex", sep="")), clean = TRUE)
  alignement <- msaConvert(multiple_alignement)
  d <-  as.matrix(dist.alignment(alignement, matrix = "identity" ))
  tree <- bionj(d)
  df <- data.frame(taxa = top_cyp_blast[i],orf = "orf")
  row.names(df) <- NULL
  ggtree(tree) %<+% df + geom_tiplab(aes(color=orf), font=2) + xlim(0,1) + geom_treescale() + theme(legend.position = "none")
  ggsave(file.path(PhylogenyDir, paste("phylo_tree_", list_cyp_phylo[i],"Rstudio.pdf", sep = "")), 
         device = "pdf",height = 60, width = 35, units = "cm", limitsize = FALSE)
}


#-------------------------------------------------------------------
#  APPENDICE : SESSION INFO                          
#-------------------------------------------------------------------
#R version 4.0.3 (2020-10-10)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19041)

#Matrix products: default

#locale:
# LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252
# LC_NUMERIC=C                   LC_TIME=French_France.1252    

#attached base packages:
# tools     stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# ape_5.5             msa_1.22.0          LncFinder_1.1.4     seqinr_4.2-5        rBLAST_0.99.2      
# Biostrings_2.58.0   XVector_0.30.0      IRanges_2.24.1      S4Vectors_0.28.1    BiocGenerics_0.36.1

# loaded via a namespace (and not attached):
# pkgload_1.2.1        splines_4.0.3        foreach_1.5.1        prodlim_2019.11.13   assertthat_0.2.1    
# BiocManager_1.30.12  remotes_2.3.0        progress_1.2.2       sessioninfo_1.1.1    ipred_0.9-11        
# pillar_1.6.0         lattice_0.20-44      glue_1.4.2           pROC_1.17.0.1        colorspace_2.0-0    
# recipes_0.1.16       Matrix_1.2-18        plyr_1.8.6           timeDate_3043.102    pkgconfig_2.0.3     
# devtools_2.4.0       caret_6.0-86         zlibbioc_1.36.0      purrr_0.3.4          scales_1.1.1        
# processx_3.5.1       gower_0.2.2          lava_1.6.9           proxy_0.4-25         tibble_3.1.1        
# generics_0.1.0       ggplot2_3.3.3        usethis_2.0.1        ellipsis_0.3.1       cachem_1.0.4        
# withr_2.4.2          nnet_7.3-15          cli_2.5.0            survival_3.2-11      magrittr_2.0.1      
# crayon_1.4.1         memoise_2.0.0        ps_1.6.0             fs_1.5.0             fansi_0.4.2         
# nlme_3.1-152         MASS_7.3-53.1        class_7.3-18         pkgbuild_1.2.0       data.table_1.14.0   
# prettyunits_1.1.1    hms_1.0.0            lifecycle_1.0.0      stringr_1.4.0        munsell_0.5.0       
# callr_3.7.0          e1071_1.7-6          ade4_1.7-16          compiler_4.0.3       rlang_0.4.10        
# grid_4.0.3           iterators_1.0.13     rstudioapi_0.13      testthat_3.0.2       ModelMetrics_1.2.2.2
# gtable_0.3.0         codetools_0.2-18     DBI_1.1.1            curl_4.3             reshape2_1.4.4      
# R6_2.5.0             lubridate_1.7.10     dplyr_1.0.5          fastmap_1.1.0        utf8_1.2.1          
# rprojroot_2.0.2      desc_1.3.0           stringi_1.5.3        Rcpp_1.0.6           vctrs_0.3.8         
# rpart_4.1-15         tidyselect_1.1.1   
