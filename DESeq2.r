library(DESeq2)
library(dplyr)

#Read in count information.
countData = read.table("Count.filtered.tsv",header=TRUE,sep = '\t')
#X an Y gene names can be the same. This makes them unique. Row is set to this unique value. 
rownames(countData)= paste(countData$Gene_name,countData$Gene_stable_ID, countData$chr,sep="$")

#Get SARs transcripts
SARS=countData[countData$chr == "SARSCOV2_ASM985889v3", ] 
#Extraxt EB and annotated genes separately. Combine with SARS
EB=countData[countData$status == "novel", ] 
ANN=countData[countData$status == "annotated", ] 
ANN_nc=ANN[ANN$Gene_type != "protein_coding", ] 
ANN_pc=ANN[ANN$Gene_type == "protein_coding", ] 

EB_Counts=rbind(EB, SARS)
ANN_Counts=rbind(ANN, SARS)
ANN_nc_Counts=rbind(ANN_nc, SARS)
ANN_pc_Counts=rbind(ANN_pc, SARS)
#Select only ID columns that comtain counts. 
ANN_Counts=select(ANN_Counts,contains("ORF10"))
EB_Counts=select(EB_Counts,contains("ORF10"))
ANN_nc_Counts=select(ANN_nc_Counts,contains("ORF10"))
ANN_pc_Counts=select(ANN_pc_Counts,contains("ORF10"))
#Round to nearest int
ANN_Counts=round(ANN_Counts,0) 
EB_Counts=round(EB_Counts,0) 
ANN_nc_Counts=round(ANN_nc_Counts,0) 
ANN_pc_Counts=round(ANN_pc_Counts,0) 
##Read in expermental design
metadata = read.table("Design.tsv",header=TRUE,row.names=1,sep = '\t')

#Should return TRUE
#all(rownames(metadata) == colnames(EB_Counts))
 




# ##############################EB
# ##Make DEseq2 object
# dds = DESeqDataSetFromMatrix(countData = EB_Counts,colData = metadata,design = ~ Batch + Condition)
# dds = DESeq(dds)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","Cont","EVC"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"ContvsEVC_EB_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","Cont","MUT"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"ContvsMUT_EB_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","Cont","WT"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"ContvsWT_EB_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","WT","MUT"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"WTvsMUT_EB_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","MUT","EVC"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"MUTvsEVC_EB_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","WT","EVC"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"WTvsEVC_EB_DGE.tsv" ,sep = '\t',row.names = FALSE)

# ##############################Annotated
# ##Make DEseq2 object
# dds = DESeqDataSetFromMatrix(countData = ANN_Counts,colData = metadata,design = ~ Batch + Condition)
# dds = DESeq(dds)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","Cont","EVC"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"ContvsEVC_Annotated_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","Cont","MUT"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"ContvsMUT_Annotated_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","Cont","WT"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"ContvsWT_Annotated_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","WT","MUT"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"WTvsMUT_Annotated_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","MUT","EVC"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"MUTvsEVC_Annotated_DGE.tsv" ,sep = '\t',row.names = FALSE)

# #Contrast case vs control
# result = results(dds, contrast=c("Condition","WT","EVC"))
# ## Remove rows with NA
# result = result[complete.cases(result),]
# #Put GeneID as column
# result = cbind(ID = rownames(result), result)
# write.table(result,"WTvsEVC_Annotated_DGE.tsv" ,sep = '\t',row.names = FALSE)

##############################Annotated  nc
##Make DEseq2 object
dds = DESeqDataSetFromMatrix(countData = ANN_nc_Counts,colData = metadata,design = ~ Batch + Condition)
dds = DESeq(dds)

#Contrast case vs control
result = results(dds, contrast=c("Condition","Cont","EVC"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"ContvsEVC_Annotated_nc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","Cont","MUT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"ContvsMUT_Annotated_nc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","Cont","WT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"ContvsWT_Annotated_nc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","WT","MUT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"WTvsMUT_Annotated_nc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","MUT","EVC"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"MUTvsEVC_Annotated_nc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","WT","EVC"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"WTvsEVC_Annotated_nc_DGE.tsv" ,sep = '\t',row.names = FALSE)

##############################Annotated pc
##Make DEseq2 object
dds = DESeqDataSetFromMatrix(countData = ANN_pc_Counts,colData = metadata,design = ~ Batch + Condition)
dds = DESeq(dds)

#Contrast case vs control
result = results(dds, contrast=c("Condition","Cont","EVC"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"ContvsEVC_Annotated_pc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","Cont","MUT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"ContvsMUT_Annotate_pcd_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","Cont","WT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"ContvsWT_Annotated_pc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","WT","MUT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"WTvsMUT_Annotated_pc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","MUT","EVC"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"MUTvsEVC_Annotated_pc_DGE.tsv" ,sep = '\t',row.names = FALSE)

#Contrast case vs control
result = results(dds, contrast=c("Condition","WT","EVC"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ID = rownames(result), result)
write.table(result,"WTvsEVC_Annotated_pc_DGE.tsv" ,sep = '\t',row.names = FALSE)

