#!/usr/bin/Rscript
# Script para generar tablas de conteo en el grupo 'case'
# APG
# Mon Sep 16 14:35:37 2019

###################################################################################################
##### Working Dir [COMENTAR/DESCOMENTAR]!!!
###################################################################################################

#workDir <- file.path("/home/amperez/Documents/TFM")
#sampleInfoDir <- file.path(paste0(workDir,"/sampleInfo"))
#dataDir <- file.path(paste0(workDir,"/rawData"))
#processDir <- file.path(paste0(workDir,"/process"))
#resultsDir <- file.path(paste0(workDir,"/results"))
#setwd(workDir)

###################################################################################################
##### Package Installation [COMENTAR/DESCOMENTAR]!!!
###################################################################################################

##¡Descomentar al lanzarlo en el cluster!
#if (!require("reshape")) install.packages("reshape") #para unir dataframes
#if (!require("data.table")) install.packages("data.table")


###################################################################################################
##### data loading [COMENTAR/DESCOMENTAR]!!!
###################################################################################################

##¡Comentar al lanzarlo en el cluster
###Crea una lista con los nombres de todos los archivos localizados en ese directorio
#results_files<-list.files("./rawData/ENOD_studies/priorizaciones_local/")
##¡Descomentar al lanzarlo en el cluster!
results_files<-list.files("./var_tables/")

###Rellenamos la lista con los archivos contenidos en los directorios
res_df<-NULL
for (i in 1:length(results_files)){
  ##¡Comentar al lanzarlo en el cluster
  #res_df[[i]]<-read.csv(file.path(paste0("./rawData/ENOD_studies/priorizaciones_local/", results_files[i])), header = TRUE, sep = "\t",stringsAsFactors=FALSE)
  ##¡Descomentar al lanzarlo en el cluster!
  res_df[[i]]<-read.csv(file.path(paste0("./var_tables/", results_files[i])), header = TRUE, sep = "\t",stringsAsFactors=FALSE)
}

###obtenemos el listado de todos los individuos
lista_ind<-NULL
for (i in 1:length(res_df)){
  ind<-as.vector(names(res_df[[i]][5]))
  lista_ind<-c(lista_ind,ind)
}
lista_ind<-sort(unique(lista_ind))
print(paste0("El número de individuos en este estudio es de: ",length(lista_ind)))

###################################################################################################
##### data modification 
###################################################################################################

###añadimos una columna con el nombre del individuo
for (i in 1:length(res_df)){
  ind_label<-rep(names(res_df[[i]][5]),length(res_df[[i]]$Variante))
  res_df[[i]]<-cbind.data.frame(res_df[[i]],ind_label)
}

###unimos todas las data.frames en 'all_df'
library(reshape)
data <- merge_recurse(res_df)
all_df<-data[,1:length(lista_ind)]

###################################################################################################
##### variant selection [COMENTAR/DESCOMENTAR]!!!
###################################################################################################

###seleccionamos variantes en funcion de CADD
selected_variants<-subset(all_df,ScaledCADD>25)

###seleccionamos variantes en regiones codificantes
coding_variants<-subset(selected_variants,Consecuencia=="coding_sequence_variant"|Consecuencia=="feature_truncation"|Consecuencia=="feature_elongation"|Consecuencia=="frameshift_variant"
                        |Consecuencia=="incomplete_terminal_codon_variant"|Consecuencia=="inframe_deletion"|Consecuencia=="inframe_insertion"|Consecuencia=="missense_variant"
                        |Consecuencia=="NMD_transcript_variant"|Consecuencia=="protein_altering_variant"|Consecuencia=="synonymous_variant"|Consecuencia=="start_lost"|Consecuencia=="stop_gained"
                        |Consecuencia=="stop_lost"|Consecuencia=="stop_retained_variant"
                        |Consecuencia=="splice_acceptor_variant"|Consecuencia=="splice_donor_variant"|Consecuencia=="transcript_ablation"|Consecuencia=="transcript_amplification")

print(paste0("Se han identificado ", nrow(coding_variants), " variantes diferentes en las regiones codificantes de los ",length(lista_ind)," individuos del estudio"))
coding_variants$Genes <- as.character(coding_variants$Genes)
####vemos qué variantes afectan a varios genes
coding_variants_dup <- coding_variants[grep(",",coding_variants$Genes),]
print(paste0("De las cuales, ", nrow(coding_variants_dup), " afectan a varios genes"))
dup <- which(coding_variants$Genes%in%coding_variants_dup$Genes)
coding_variants_uniq <- coding_variants[-dup,]
CADD_cases_uniq <- cbind.data.frame(coding_variants_uniq$Genes,round(coding_variants_uniq$ScaledCADD,2))
colnames(CADD_cases_uniq) <- c("genes","CADD")
CADD_cases_dup <- cbind.data.frame(coding_variants_dup$Genes,round(coding_variants_dup$ScaledCADD,2))
colnames(CADD_cases_dup) <- c("genes","CADD") 

#[COMENTAR/DESCOMENTAR]!!!
#añadir dir cluster
#write.csv(CADD_cases_uniq, paste0(resultsDir,"/CADD_cases_uniq.csv"),row.names = F)
#write.csv(CADD_cases_dup, paste0(resultsDir,"/CADD_cases_dup.csv"),row.names = F)
write.csv(CADD_cases_uniq, "./results/CADD_cases_uniq.csv",row.names = F)
write.csv(CADD_cases_dup, "./results/CADD_cases_dup.csv",row.names = F)

###cuántas variantes hay en cada gen?
VARperGENE_dup <- table(coding_variants_dup$Genes)
VARperGENE_dup <- as.data.frame(VARperGENE_dup)
VARperGENE_uniq <- table(coding_variants_uniq$Genes)
VARperGENE_uniq <- as.data.frame(VARperGENE_uniq)
#write.csv(VARperGENE_dup, paste0(resultsDir,"/VARperGENE_dup.csv"),row.names = F)
#write.csv(VARperGENE_uniq, paste0(resultsDir,"/VARperGENE_uniq.csv"),row.names = F)
write.csv(VARperGENE_dup, "./results/VARperGENE_dup.csv",row.names = F)
write.csv(VARperGENE_uniq, "./results/VARperGENE_uniq.csv",row.names = F)

###################################################################################################
##### Conteos [COMENTAR/DESCOMENTAR]!!!
###################################################################################################

###obtenemos el listado de genes que contienen las variantes seleccionadas
lista_genes_uniq<-sort(unique(coding_variants$Genes))
lista_genes_uniq<-lista_genes_uniq[!lista_genes_uniq==""]
#write.table(lista_genes_uniq,paste0(resultsDir,"/lista_genes_uniq.txt"),row.names = F, col.names = F, quote = F)
write.table(lista_genes_uniq,"./results/lista_genes_uniq.txt",row.names = F, col.names = F, quote = F)
#lg<-paste(lista_genes,collapse=",")
#lista_genes2<-unlist(strsplit(lg, split=","))
lista_genes_all<-sort(coding_variants$Genes)
lista_genes_all<-lista_genes_all[!lista_genes_all==""]
#write.table(lista_genes_all,paste0(resultsDir,"/lista_genes_all.txt"),row.names = F, col.names = F, quote = F)
write.table(lista_genes_all,"./results/lista_genes_all.txt",row.names = F, col.names = F, quote = F)
  
###separamos en tantas data.frames como ind haya
filtered_df<-lapply(unique(coding_variants$ind_label), function(x) coding_variants[coding_variants$ind_label == x,])

###listado de ind
#comprobamos que el orden de los ind es el mismo en la tabla conjunta que en los df por separado
listado <- NULL #en la tabla conjunta
for (i in 1:length(lista_ind)){
  ind<- unique(as.character(coding_variants$ind_label))[i]
  listado <- c(listado,ind)
}
listado2 <- NULL #en los data frames por separado
for (i in 1:length(filtered_df)){
  ind2<- as.character(unique(filtered_df[[i]]$ind_label))
  listado2 <- c(listado2,ind2)
}
print(paste0("El orden de los individuos en el df conjunto es igual que en la lista con los df por separado: ",all(listado==listado2)))

###funcion para obtener listado de genes por ind
genes_individuos <- function(genes){
  gen_listado <- genes$Genes
  gen_listado <- paste(gen_listado,collapse=",")
  gen_listado <- unlist(strsplit(gen_listado, split=","))
  return(gen_listado)
}
#ejecutamos la funcion en el objeto con lo df pro separado
gen_lists <- lapply(filtered_df,genes_individuos)
#nombramos cada uno de los listados con el label del individuo
names(gen_lists) <- c(listado)

for(i in 1:length(gen_lists)){
  #cambiar directorio para cluster!!!
  write.table(gen_lists[[i]],file=file.path(paste0("./results/",names(gen_lists)[[i]],"_genes_label.txt")),row.names = F,col.names = names(gen_lists)[[i]])
  write.table(gen_lists[[i]],file=file.path(paste0("./results/",names(gen_lists)[[i]],"_genes.txt")),row.names = F,col.names = F)
  
}

for(i in 1:length(gen_lists)){
  write.csv(gen_lists[[i]],file=file.path(paste0("./results/",names(gen_lists)[[i]],"_genes.csv")),row.names = F)
}

suma<-NULL
for(i in 1:length(gen_lists)){
  #cambiar directorio para cluster!!!
  size<-length(gen_lists[[i]])
  suma<-sum(size,suma)
}

print(paste0("El total de genes entre todos los individuos es de: ",suma))


#creamos la matriz indxgenes vacía
bin_cases<-matrix(nrow = 0, ncol = length(lista_genes_uniq))
colnames(bin_cases)<-c(lista_genes_uniq)

#la rellenamos consultando en la lista de genes de cada individuo 
for(i in 1:length(filtered_df)){
  g<-unique(filtered_df[[i]]$Genes)
  gvar<-as.numeric(lista_genes_uniq%in%g)
  bin_cases<-rbind(bin_cases,gvar)
}
rownames(bin_cases)<-c(lista_ind)

bin_cases_dup <- bin_cases
bin_cases_dup <- bin_cases_dup[,grep(",",colnames(bin_cases_dup))]
rownames(bin_cases_dup)<-c(lista_ind)
dup <- which(colnames(bin_cases)%in%colnames(bin_cases_dup))
bin_cases_uniq <- bin_cases[,-dup]
rownames(bin_cases_uniq)<-c(lista_ind)

#write.csv(bin_cases_dup, paste0(resultsDir,"/bin_cases_dup.csv"))
#write.csv(bin_cases_uniq, paste0(resultsDir,"/bin_cases_uniq.csv"))
write.csv(bin_cases_dup, "./results/bin_cases_dup.csv")
write.csv(bin_cases_uniq, "./results/bin_cases_uniq.csv")



#generamos la tabla de contingencia de los casos
counts_cases_dup<-apply(bin_cases_dup,2,table)
#write.csv(counts_cases_dup,paste0(resultsDir,"/counts_cases_dup.csv"))
counts_cases_uniq<-apply(bin_cases_uniq,2,table)
#write.csv(counts_cases_uniq,paste0(resultsDir,"/counts_cases_uniq.csv"))
write.csv(counts_cases_dup,"./results/counts_cases_dup.csv")
write.csv(counts_cases_uniq,"./results/counts_cases_uniq.csv")

##guardamos todo
save.image("cases.RData")
