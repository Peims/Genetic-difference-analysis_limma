rm(list=ls())
a<- data.table::fread('miRNA.txt_1') #%>% dplyr::select(-c(2:4))
library(tibble)
a<- column_to_rownames(a,'sRNA.readcount')

library(limma)

class <- c(rep("con",3),rep("treat",3))
design <- model.matrix(~factor(class))
colnames(design) <- c("con","treat")
design
#线性模型拟合
fit <- lmFit(a,design)
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
write.table(allDiff,'miRNAdiff_out',sep="\t",row.names = T,quote=F)
