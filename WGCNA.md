# Building co-expression network using WGCNA
BioProject PRJEB25639 contained 69 samples with three biological replicates and one sample with two biological replicates.
We use these data to constructed co-expression network.
Firstly, we got gene level counts and then got normolized counts using DESeq2.
####Got normolized counts and log2(n+1) transformation
```R
# working dir is in PRJEB25639
library("DESeq2")
library(tximportData)
library(tximport)
library(readr)
samples <- read.table(file.path("./salmon_quant", "sample_information.txt"), header = TRUE, sep = '\t')
samples
files <- file.path("./salmon_quant", samples$File, "quant.sf")
files
names(files) <- paste0(samples$sample)
samples$condition <- factor(samples$Order)
all(file.exists(files))
t2g <- read.table(file.path("/data2/user_data/rna_seq", "gene_transcript_list.txt"), header = TRUE)
t2g <- dplyr::rename(t2g,target_id = transcript_id, gene_id = gene_id)
head(t2g)
txi <- tximport(files, type = "salmon", tx2gene = t2g)
dds <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
nor_counts <- counts(dds, normalized=TRUE)
write.table(as.data.frame(nor_counts), file="./salmon_quant/BCS_normalized_count.tsv")
# this gives log2(n + 1)
ntd <- normTransform(dds)
write.table(as.data.frame(assay(ntd)), file="./salmon_quant/BCS_normalized_count_log.tsv")
```
####genes with tpm ≥ 5 at least 10 samples, and genes with tpm < 0.5 was view as not expressed.


```R
dirpath="/data2/masw_data/transcript/TSG/WGCNA/"
setwd(dirpath) 
library(WGCNA)
enableWGCNAThreads(nThreads = 30)
set.seed(1)
options(stringsAsFactors = F)
input_origin<-read.csv("WGCNA_input_more_than_5FPKM.txt",head=TRUE,row.names="Geneid",sep="\t")
input=t(input_origin)
nGenes=ncol(input)
nSamples=nrow(input)
datExpre_tree<-hclust(dist(input),method="average")


#########################查看sample聚类情况####################################
pdf("Sample_cluster2.pdf",width=20,height=20)
plot(datExpre_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 4, 
     cex.axis = 2, cex.main = 4,cex.lab=2)
dev.off()

###################################着色#########################################
sample_colors <- numbers2colors(input, colors = c("white","blue"),signed = FALSE)
pdf("Sample_heat_cluster.pdf",width=25,height=60)
plotDendroAndColors(datExpre_tree, sample_colors,
                    groupLabels = colnames(input),
                    cex.dendroLabels = 0.8,#聚类的sample的文字
                    marAll = c(4, 16, 12,4),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")
dev.off()
####################################主元分析####################################
png("PCA.png")
pca = prcomp(input)
sampletype<-rownames(input)
plot(pca$x[,c(1,2)],pch=16,col=rep(rainbow(nSamples),each=1),cex=1.5,main = "PCA map")
text(pca$x[,c(1,2)],row.names(pca$x),col="black",pos=3,cex=1)
legend("right",legend=sampletype,ncol = 1,xpd=T,inset = -0.15,
       pch=16,cex=1,col=rainbow(length(sampletype)),bty="n")
dev.off()
###########################计算软阈值####################################################
powers=c(1:30)
pow<-pickSoftThreshold(input, powerVector = powers, verbose = 5)

#####################设置网络构建参数选择范围，计算无尺度分布拓扑矩阵#####################
pdf("softpower.pdf",width=10,height=10)
plot(pow$fitIndices[,1], -sign(pow$fitIndices[,3])*pow$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(pow$fitIndices[,1], -sign(pow$fitIndices[,3])*pow$fitIndices[,2],labels=powers,cex=0.5,col="red")
plot(pow$fitIndices[,1], pow$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(pow$fitIndices[,1], pow$fitIndices[,5], labels=powers, cex=0.6,col="red")
dev.off()
#####################构架基因网络###########################################################
  write.csv(pow,"softPower.select.csv")
  softPower=11
  dir.create(paste(dirpath,softPower,sep=""))
  adjacency = adjacency(input, power = softPower) ##计算树之间的邻接性
  TOM = TOMsimilarity(adjacency)##计算树之间的相异性
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average") ##聚类分析
  
  minModuleSize = 30  #设置基因模块中最少基因包含的数目
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = TRUE,pamRespectsDendro = FALSE, minClusterSize = minModuleSize)  #分类
  
  dynamicColors = labels2colors(dynamicMods)  ##分组上色
  table(dynamicColors)
  
  
  ##################计算基因模块的特征值###############################
  MEList = moduleEigengenes(input, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs)
  METree = hclust(as.dist(MEDiss), method = "average")
  
  #########################建立系统聚类树#################################
  MEDissThres = 0.1  #可改，类似于bootstrap值的东西,越低，聚类的要求就越高
  pdf(paste(softPower,"/Gene.cluster.pdf",sep=""),width=10,height=10)
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  
  ######################基因模块合并###########################
  
  merge = mergeCloseModules(input, dynamicColors, MEs = MEs,cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors
  mergedMEs = merge$newMEs
  write.csv(MEDiss,paste(softPower,"/MEDiss.csv",sep=""))
  write.csv(MEs,paste(softPower,"/MEs.csv",sep=""))
  
  #########################绘制基因网络图###########################################
  pdf(paste(softPower,"/geneTree.pdf",sep=""),width=10,height=10)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  
  #########################不同模块基因热图及关键基因的表达############################
  person=cor(input,use = 'p')
  
  corr<-TOM
  colnames(corr)<-colnames(input)
  rownames(corr)<-colnames(input)
  Colors<-mergedColors
  names(Colors)<-colnames(input)
  
  colnames(person)<-colnames(input)
  rownames(person)<-colnames(input)
  
  write.csv(Colors,paste(softPower,"/Colors.csv",sep=""))
  
  umc = unique(mergedColors)
  lumc = length(umc)
  
  for (i in c(1:lumc)){
    if(umc[i]== "grey"){
      next
    }
    ME=MEs[, paste("ME",umc[i], sep="")]
    outfile1=paste(softPower,"/Diff_module.",i,".heatmap.png",sep="")
    outfile2=paste(softPower,"/Diff_module.",i,".expression.png",sep="");
    png(outfile1)
    plotMat(t(scale(input[,Colors==umc[i]])),nrgcols=30,rlabels=F,rcols=umc[i], main=umc[i], cex.main=2)
    dev.off()
    png(outfile2)
    barplot(ME, col=umc[i], main="", cex.main=2,ylab="eigengene expression",xlab="array sample")
    dev.off()
  }
  ##############################基因共表达网络热图#######################
  kME=signedKME(input, mergedMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
  # 简称MM, 将该基因的表达量与module的第一主成分，即module eigengene进行相关性分析就可以得到MM值，所以MM值本质上是一个相关系数，如果基因和某个module的MM值为0，说明二者根本不相关，该基因不属于这个module; 如果MM的绝对值接近1，说明基因与该module相关性很高。
  head(kME)
  if (dim(input)[2]>=1500) nSelect=1500 else nSelect=dim(input)[2]
  set.seed(1)
  select = sample(nGenes, size = nSelect)
  selectTOM = dissTOM[select, select]
  selectTree = hclust(as.dist(selectTOM), method = "average")
  mergedColors = merge$colors
  moduleColors = mergedColors  #origin :"moduleColors=lables2colors(mergedColors)" fixed, no meaning except for make myself misunderstanding
  selectColors = moduleColors[select]
  pdf(paste(softPower,"/network.heatmap.plot.pdf",sep=""),width=10,height=10)
  plotDiss = selectTOM^7
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot")
  dev.off()
  ##################导出网络到Cytoscape#######################################
  #选择需要导出的模块颜色
  for(modules in umc){
    probes = colnames(input)
    inModule = is.finite(match(moduleColors, modules));
    modProbes = probes[inModule];
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    # Export the network into edge and node list files Cytoscape can read
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = paste(softPower,"/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                   nodeFile = paste(softPower,"/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = 0.02,
                                   nodeNames = modProbes,
                                   #altNodeNames = modGenes,
                                   nodeAttr = moduleColors[inModule]);
    
  }

  ```