library(NMF)
library(dplyr)
library(limma)
library(R.matlab)
source("/your_path_to/io.R")
### control information
untrt.info <- read.delim2("/your_path_to/inst.info",header = T,stringsAsFactors = F,check.names = F) %>%
  dplyr::select(.,one_of(c("distil_id","pert_id","pert_desc","cell_id","pert_time","pert_dose","pert_dose_unit",
                           "pert_type"))) %>%
  tbl_df %>%
  filter(.,pert_type=="ctl_vehicle")%>%
  filter(.,cell_id=="HEPG2")%>%
  filter(.,pert_desc!="-666")
untrt.infomat <- as.matrix(untrt.info)
nr <- nrow(untrt.infomat)
nc <- ncol(untrt.infomat)
untrt.infomat1 <- matrix(0,nr,(nc+1))
untrt.infomat1[,1:nc] <- untrt.infomat
for (i in 1:nr)
{
  A <- strsplit(untrt.infomat1[i,1],":")
  untrt.infomat1[i,(nc+1)] <- A[[1]][1]
}
message(1)

### read the GPL96 annotation data
GPL96_anno <- read.delim2("/your_path_to/GPL961.annot",header = T)
GPL96_anno <- as.matrix(GPL96_anno)
### matching probe data to gene symbol data
## matching probe ids to gene symbols
probe_to_gene_96 <- matrix(0,54675,10)
probes_96 <- parse.gctx("/your_path_to/q2norm_n1328098x22268.gctx",cid=untrt.infomat1[1,1])@rid
gene_symbol_96 <- vector("character")
count1 <- 0
for (i in 1:length(probes_96))
{
  B <- GPL96_anno[which(GPL96_anno[,1]==probes_96[i],T),"Gene_symbol"]
  if (B!="")
  {
    probe_to_gene_96[i,1] <- probes_96[i]
    A <- strsplit(B,"///")[[1]]
    probe_to_gene_96[i,2:(length(A)+1)] <- A
    count1 <- count1+1
    gene_symbol_96 <- c(gene_symbol_96,A)
  }
}
probe_to_gene_96 <- probe_to_gene_96[1:count1,]
gene_symbol_96 <- unique(gene_symbol_96)
gene_to_probe_96 <- matrix(0,length(gene_symbol_96),10)
for (i in 1:length(gene_symbol_96))
{
  gene_to_probe_96[i,1] <- gene_symbol_96[i]
  matched_probes_96 <- probe_to_gene_96[which(probe_to_gene_96==gene_symbol_96[i],T)[,1],1]
  gene_to_probe_96[i,2:(length(matched_probes_96)+1)] <- matched_probes_96
}
message(2)

### chemical compound treated information
inst.info <- read.delim2("/your_path_to/inst.info",header = T,stringsAsFactors = F,check.names = F) %>%
  dplyr::select(.,one_of(c("distil_id","pert_id","pert_desc","cell_id","pert_time","pert_dose","pert_dose_unit",
                           "pert_type"))) %>%
  tbl_df %>%
  filter(.,pert_type=="trt_cp")%>%
  filter(.,cell_id=="HEPG2")%>%
  filter(.,pert_desc!="-666")
inst.infomat <- as.matrix(inst.info)
nr <- nrow(inst.infomat)
nc <- ncol(inst.infomat)
inst.infomat1 <- matrix(0,nr,(nc+1))
inst.infomat1[,1:nc] <- inst.infomat
for (i in 1:nr)
{
  A <- strsplit(inst.infomat1[i,1],":")
  inst.infomat1[i,(nc+1)] <- A[[1]][1]
}
message(3)

### classify info by the combination of cp, dose and time
unique_info <- unique(inst.infomat[,c("pert_id","pert_time","pert_dose")])
info_list <- vector("list",nrow(unique_info))
for (i in 1:length(info_list))
{
  sub_list <- vector("list",5)
  names(sub_list) <- c("pert_id","pert_time","pert_dose","treated_distil_id","control_distil_id")
  sub_list[[1]] <- unique_info[i,1]
  sub_list[[2]] <- unique_info[i,2]
  sub_list[[3]] <- unique_info[i,3]
  sub_info <- inst.infomat1[which(inst.infomat1[,"pert_id"]==unique_info[i,1],T),]
  sub_info <- sub_info[which(sub_info[,"pert_time"]==unique_info[i,2],T),]
  sub_info <- sub_info[which(sub_info[,"pert_dose"]==unique_info[i,3],T),]
  sub_list[[4]] <- sub_info[,1]
  batch_info <- unique(sub_info[,(nc+1)])
  ctl_info <- untrt.infomat1[which(untrt.infomat1[,(nc+1)]%in%batch_info,T),]
  sub_list[[5]] <- ctl_info[,1]
  info_list[[i]] <- sub_list
}
message(4)

### extract data to calculate fold change with limma
lincs_fold_change <- matrix(0,22268,length(info_list))
lincs_probe_id <- parse.gctx("/your_path_to/q2norm_n1328098x22268.gctx",cid=inst.infomat1[1,1])@rid
lincs_col_name <- vector("character",length(info_list))
for (i in 1:length(info_list))
{
  ## extract cp treated data
  sample_info <- info_list[[i]][[4]]
  sample_data <- matrix(0,22268,length(sample_info))
  for (j in 1:ncol(sample_data))
  {
    sample_data[,j] <- parse.gctx("/your_path_to/q2norm_n1328098x22268.gctx",cid=sample_info[j])@mat
  }
  ## extract control data
  control_info <- info_list[[i]][[5]]
  control_data <- matrix(0,22268,length(control_info))
  for (j in 1:ncol(control_data))
  {
    control_data[,j] <- parse.gctx("/your_path_to/q2norm_n1328098x22268.gctx",cid=control_info[j])@mat
  }
  ## transfer the probe data to gene symbol data
  new_data1 <- cbind(sample_data,control_data)
  new_data <- matrix(0,length(gene_symbol_96),ncol(new_data1))
  for (i in 1:length(gene_symbol_96))
  {
    sub_probes <- gene_to_probe_96[i,which(gene_to_probe_96[i,-1]!="0",T)]
    if (length(sub_probes)>1)
    {
      sub_matr1 <- new_data1[which(probes_96%in%sub_probes,T),]
      new_sub_matr1 <- apply(sub_matr1,2,as.numeric)
      new_data[i,] <- apply(new_sub_matr1,2,mean)
    }
    else
    {
      sub_matr1 <- new_data1[which(probes_96%in%sub_probes,T),]
      new_data[i,] <- apply(sub_matr1,as.numeric)
    }
  }
  rownames(new_data) <- gene_symbol_96
  ## t test with limma
  pheno = data.frame(state=c(rep("cancer",ncol(sample_data)),rep("normal",ncol(control_data))))
  pheno$state = as.factor(pheno$state)
  design <- model.matrix(~ pheno$state + 0)
  colnames(design) = levels(pheno$state)
  fit <- lmFit(new_data, design)
  colnames(design) = levels(pheno$state)
  cont.matrix <- makeContrasts(cancer-normal, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  tT <- topTable(fit2, adjust="BH",sort.by="B",number=dim(new_data)[1])
  tT$ID <- rownames(tT)
  ## assign the t value as fold change 
  t_value <- vector("numeric",nrow(new_data))
  for (j in 1:nrow(tT))
  {
    index_j <- which(tT$ID==gene_symbol_96[j],T)
    t_value[j] <- tT$t[index_j]
  }
  lincs_fold_change[,i] <- t_value
  lincs_col_name[i] <- paste(info_list[[i]][1],"---",info_list[[i]][2],"---",info_list[[i]][3],"---")
}
rownames(lincs_fold_change) <- unlist(gene_symbol_96)
colnames(lincs_fold_change) <- unlist(lincs_col_name)
write.table(gene_symbol_96,"rownames_of_lincs_fold_change_data_by_gene_symbol.txt",sep="\t",quote=T,row.names = F,col.names = F)
write.table(lincs_col_name,"colnames_of_lincs_fold_change_data_by_gene_symbol.txt",sep="\t",quote=T,row.names = F,col.names = F)
writeMat("lincs_fold_change_data_by_gene_symbol.mat",A=lincs_fold_change)
message(5)

### read the geo data and classify by experiment time
c3r1 <- read.delim2("/your_path_to/GSE72068_series_matrix_processed.txt",header = T)
c3r1 <- as.matrix(c3r1)
c3_data <- c3r1[,c(9,10,11,20)]
rownames(c3_data) <- c3r1[,1]
write.table(c3_data,"geo_hepatitis_data_with_12_days_by_probe.txt",sep="\t",quote=F,row.names=T,col.names = T)
message(6)

### read in the annotation data of GPL570 plateform
anno <- read.delim2("/your_path_to/GPL10558-50081.txt",header=T)
anno <- as.matrix(anno)
message(7)

### matching probe data to gene symbol data
## matching probe ids to gene symbols
probe_to_gene_570 <- matrix(0,54675,10)
probes_570 <- c3r1[,1]
gene_symbol_570 <- vector("character")
count1 <- 0
for (i in 1:length(probes_570))
{
  B <- anno[which(anno[,1]==probes_570[i],T),"Gene.Symbol"]
  if (B!="")
  {
    probe_to_gene_570[i,1] <- probes_570[i]
    A <- strsplit(B," /// ")[[1]]
    probe_to_gene_570[i,2:(length(A)+1)] <- A
    count1 <- count1+1
    gene_symbol_570 <- c(gene_symbol_570,A)
  }
}
probe_to_gene_570 <- probe_to_gene_570[1:count1,]
gene_symbol_570 <- unique(gene_symbol_570)
gene_to_probe_570 <- matrix(0,length(gene_symbol_570),10)
for (i in 1:length(gene_symbol_570))
{
  gene_to_probe_570[i,1] <- gene_symbol_570[i]
  matched_probes_570 <- probe_to_gene_570[which(probe_to_gene_570==gene_symbol_570[i],T)[,1],1]
  gene_to_probe_570[i,2:(length(matched_probes_570)+1)] <- matched_probes_570
}
## matching probe data to gene symbol data
c3_gene_data <- matrix(0,length(gene_symbol_570),ncol(c3_data))
for (i in 1:length(gene_symbol_570))
{
  sub_probes <- gene_to_probe_570[i,which(gene_to_probe_570[i,-1]!="0",T)]
  if (length(sub_probes)>1)
  {
    sub_matr1 <- c3_data[which(probes_570%in%sub_probes,T),]
    new_sub_matr1 <- apply(sub_matr1,2,as.numeric)
    c3_gene_data[i,] <- apply(new_sub_matr1,2,mean)
  }
  else
  {
    sub_matr1 <- c3_data[which(probes_570%in%sub_probes,T),]
    c3_gene_data[i,] <- apply(sub_matr1,as.numeric)
  }
}
rownames(c3_gene_data) <- gene_symbol_570
colnames(c3_gene_data) <- colnames(c3_data)
write.table(c3_gene_data,"geo_hepatitis_data_with_12_days_by_gene_symbol.txt",sep="\t",quote=T,row.names=T,col.names = T)
write.table(gene_symbol_570,"gene_symbols_in_GPL10558_platform.txt",sep="\t",quote = T,row.names = F,col.names = F)
message(8)

### t test of geo data
## t test of c3 data
pheno = data.frame(state=c(rep("cancer",2),rep("normal",2)))
pheno$state = as.factor(pheno$state)
design <- model.matrix(~ pheno$state + 0)
colnames(design) = levels(pheno$state)
fit <- lmFit(c3_gene_data, design)
colnames(design) = levels(pheno$state)
cont.matrix <- makeContrasts(cancer-normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="BH",sort.by="B",number=dim(c3_gene_data)[1])
tT$ID <- rownames(tT)
c3_tT <- as.matrix(tT)
write.table(c3_tT,"t_test_result_of_geo_hepatitis_data.txt",quote=F,sep="\t")
message(9)

### select the shared gene symbols
shared_gene_symbol <- intersect(gene_symbol_96,gene_symbol_570)
new_lincs_fold_change <- lincs_fold_change[which(gene_symbol_96%in%shared_gene_symbol,T),]
new_c3_fold_change <- c3_tT[which(c3_tT[,7]%in%shared_gene_symbol,T),3]
rownames(new_c3_fold_change) <- shared_gene_symbol
message(10)

data_normalize <- function(X)
{
  mean_X <- mean(X)
  var_X <- sqrt(var(X))
  res <- (X-mean_X)/var_X
  return(res)
}
len1 <- length(new_c3_fold_change)
all_data <- cbind(new_lincs_fold_change,as.matrix(new_c3_fold_change,len1,1))
all_data <- apply(all_data,2,data_normalize)
min_value <- min(all_data)
all_data <- all_data+abs(min_value)+10
estim.r <- nmf(all_data,2:30,nrun=10,seed=123456)


pdf(file = "factorization_rank_estimation.pdf")
plot(estim.r,width=20,height=2)
dev.off()
plot(estim.r)
message(11)

k=17
nmf_res <- nmf(t(all_data),k,nrun=10)

row_dim <- nrow(nmf_res@fit@W)
dist_matrix <- matrix(0,row_dim,row_dim)
for (i in 1:row_dim)
{
  for (j in 1:row_dim)
  {
    dist_matrix[i,j] <- sqrt(sum((nmf_res@fit@W[i,]-nmf_res@fit@W[j,])^2))
  }
}
centers <- matrix(0,k,k)
centers[1,] <- nmf_res@fit@W[row_dim,]
index <- row_dim#the index of hbv
max_index <- which(dist_matrix[index,]==max(dist_matrix[index,]),T)
centers[2,] <- nmf_res@fit@W[max_index,]
index <- c(index,max_index)
for (i in 3:k)
{
  dist11 <- apply(dist_matrix[index,],2,mean)
  max_index <- which(dist11==max(dist11),T)
  while (max_index%in%index)
  {
    dist11[max_index] <- 0
    max_index <- which(dist11==max(dist11),T)
  }
  centers[i,] <- nmf_res@fit@W[max_index,]
  index <- c(index,max_index)
}

yyy <- kmeans(nmf_res@fit@W,centers = centers,iter.max = 1000,nstart = 20)
clu <- yyy$cluster
