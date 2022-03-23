library(DESeq2)
library(dplyr)

#Read in count information.
countData = read.table("Count.filtered.tsv",header=TRUE,sep = '\t')
#X an Y gene names can be the same. This makes them unique. Row is set to this unique value. 
rownames(countData)= paste(countData$Gene_stable_ID, countData$chr,sep="_")

#Get SARs transcripts
SARS=countData[countData$chr == "SARSCOV2_ASM985889v3", ] 
#Extraxt EB and annotated genes separately
EB=countData[countData$status == "novel", ] 
ANN=countData[countData$status == "annotated", ] 


#Select only ID columns that comtain counts. 
countData=select(countData,contains("SRR"))
#Round to nearest int
countData=round(countData,0)
##Read in expermental design
metadata = read.table("Design.tsv",header=TRUE,row.names=1,sep = '\t')


#Get sample names
samples=row.names(metadata)


#Should return TRUE
#all(rownames(metadata) == colnames(countData))

##Make DEseq2 object
dds = DESeqDataSetFromMatrix(countData = countData,colData = metadata,design = ~ Sample)
dds = DESeq(dds)
#Contrast case vs control
result = results(dds, contrast=c("Sample","EVC","WT"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(miRNA_ID = rownames(result), result)


write.table(result,"EVCvsWT_DGE.tsv" ,sep = '\t',row.names = FALSE)
