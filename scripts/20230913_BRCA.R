
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RisktypeProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    scale_fill_manual(values = group_cols)+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 将图表标题居中
  return(p)
}

#筛选编码蛋白基因#########
genecode=read.delim('GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)


###TCGA#############
#######表达谱
tcga.tumor.fpkm=read.delim('origin_datas/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt',check.names = F)
tcga.tumor.fpkm[1:5,1:5]
tcga.tumor.fpkm=tcga.tumor.fpkm[which(tcga.tumor.fpkm[, 1]!=''),] 
tcga.tumor.fpkm=tcga.tumor.fpkm[!duplicated(tcga.tumor.fpkm$Hugo_Symbol),]
rownames(tcga.tumor.fpkm)=tcga.tumor.fpkm$Hugo_Symbol
tcga.tumor.fpkm=tcga.tumor.fpkm[,-c(1,2)]

tcga.normal.fpkm=read.delim('origin_datas/brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt',check.names = F)
tcga.normal.fpkm[1:5,1:5]
tcga.normal.fpkm=tcga.normal.fpkm[which(tcga.normal.fpkm[, 1]!=''),] 
tcga.normal.fpkm=tcga.normal.fpkm[!duplicated(tcga.normal.fpkm$Hugo_Symbol),]
rownames(tcga.normal.fpkm)=tcga.normal.fpkm$Hugo_Symbol
tcga.normal.fpkm=tcga.normal.fpkm[,-c(1,2)]

com.genes=intersect(rownames(tcga.tumor.fpkm),rownames(tcga.normal.fpkm))
length(com.genes)
#20486

tcga.fpkm=cbind(tcga.tumor.fpkm[com.genes,],tcga.normal.fpkm[com.genes,])
dim(tcga.fpkm)
range(tcga.fpkm)
tcga.fpkm=log2(tcga.fpkm+1)


tcga.type=data.frame(Samples=c(colnames(tcga.tumor.fpkm),colnames(tcga.normal.fpkm)),
                     Type=rep(c('Tumor','Normal'),c(length(colnames(tcga.tumor.fpkm)),length(colnames(tcga.normal.fpkm)))))
rownames(tcga.type)=tcga.type$Samples
table(tcga.type$Type)

#####临床信息
tcga.cli=read.delim('origin_datas/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018_clinical_data.tsv',check.names = F)
head(tcga.cli)
table(tcga.cli$`Cancer Type Detailed`)
tcga.cli=data.frame(Samples=tcga.cli$`Sample ID`,Age=tcga.cli$`Diagnosis Age`,
                    Stage=tcga.cli$`Neoplasm Disease Stage American Joint Committee on Cancer Code`,
                    T.stage=tcga.cli$`American Joint Committee on Cancer Tumor Stage Code`,
                    N.stage=tcga.cli$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`,
                    M.stage=tcga.cli$`American Joint Committee on Cancer Metastasis Stage Code`,
                    OS=str_split_fixed(tcga.cli$`Overall Survival Status`,':',2)[,1],
                    OS.time=as.numeric(tcga.cli$`Overall Survival (Months)`/12*365),
                    DSS=str_split_fixed(tcga.cli$`Disease-specific Survival status`,':',2)[,1],
                    DSS.time=as.numeric(tcga.cli$`Months of disease-specific survival`/12*365),
                    PFS=str_split_fixed(tcga.cli$`Progression Free Status`,':',2)[,1],
                    PFS.time=as.numeric(tcga.cli$`Progress Free Survival (Months)`/12*365))
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli[tcga.cli$OS.time>30 & tcga.cli$OS.time<3650,]

table(tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage=='STAGE X']=NA
tcga.cli$Stage=gsub('STAGE ','',tcga.cli$Stage)
tcga.cli$Stage=gsub('[ABC]','',tcga.cli$Stage)

table(tcga.cli$T.stage)
tcga.cli$T.stage[tcga.cli$T.stage=='TX']=NA
tcga.cli$T.stage=gsub('[ABCD]','',tcga.cli$T.stage)

table(tcga.cli$N.stage)
tcga.cli$N.stage[tcga.cli$N.stage=='NX']=NA
tcga.cli$N.stage[tcga.cli$N.stage=='N0 (MOL+)']='N0'
tcga.cli$N.stage[tcga.cli$N.stage=='N1MI']='N1'
tcga.cli$N.stage=gsub('[ABCI()+-]','',tcga.cli$N.stage)
tcga.cli$N.stage=gsub(' ','',tcga.cli$N.stage)

table(tcga.cli$M.stage)
tcga.cli$M.stage[tcga.cli$M.stage=='MX']=NA
tcga.cli$M.stage=gsub('[ABCI()+-]','',tcga.cli$M.stage)
tcga.cli$M.stage=gsub(' ','',tcga.cli$M.stage)

######994个样本
tcga.sample=intersect(tcga.cli$Samples,colnames(tcga.fpkm))
length(tcga.sample)
#994
range(tcga.fpkm)
tcga.exp=tcga.fpkm[intersect(rownames(tcga.fpkm),mrna_genecode$SYMBOL),tcga.sample]
dim(tcga.exp)
tcga.cli=tcga.cli[tcga.sample,]
dim(tcga.cli)

# table(tcga.cli$OS.time>30,tcga.cli$OS.time<3650)

#########GSE20685#####
load('origin_datas/GEO/GSE20685.RData')
GSE20685.cli=pData(GSE20685)
head(GSE20685.cli)
GSE20685.cli=data.frame(Samples=GSE20685.cli$geo_accession,
                          Age=GSE20685.cli$`age at diagnosis:ch1`,
                          T.stage=GSE20685.cli$`t_stage:ch1`,
                          N.stage=GSE20685.cli$`n_stage:ch1`,
                          M.stage=GSE20685.cli$`m_stage:ch1`,
                          OS.time=GSE20685.cli$`follow_up_duration (years):ch1`,
                          OS=GSE20685.cli$`event_death:ch1`)
rownames(GSE20685.cli)=GSE20685.cli$Samples
head(GSE20685.cli)
GSE20685.cli=crbind2DataFrame(GSE20685.cli)
GSE20685.cli$OS.time=GSE20685.cli$OS.time*365
GSE20685.cli$T.stage=paste0('T',GSE20685.cli$T.stage)
GSE20685.cli$N.stage=paste0('N',GSE20685.cli$N.stage)
GSE20685.cli$M.stage=paste0('M',GSE20685.cli$M.stage)

GSE20685.df=exprs(GSE20685)
GSE20685.df[1:5,1:5]
range(GSE20685.df)
GSE20685.exp=exp_probe2symbol_v2(datExpr = GSE20685.df,GPL = 'GPL570')
dim(GSE20685.exp)


######01.糖酵解在BRCA中作用########
dir.create('results/01.prognosis')
###########1.1 糖酵解与预后#####################
load('results/hall.ssGSEA.all.RData') 
my_violin(dat = tcga.hall.ssGSEA.all['HALLMARK_GLYCOLYSIS',tcga.type$Samples],group = tcga.type$Type,
          ylab = 'HALLMARK_GLYCOLYSIS',legend.position = 'none',xlab = 'Tissue',group_cols = pal_nejm()(8)[c(3,6)])+
  ggtitle('TCGA')
ggsave('results/01.prognosis/TCGA_glycolysis_boxplot.pdf',height = 5,width = 5)

load('results/hall.ssGSEA.RData')
tcga.hall.ssGSEA=tcga.hall.ssGSEA[,tcga.cli$Samples]
GSE20685.hall.ssGSEA=GSE20685.hall.ssGSEA[,GSE20685.cli$Samples]

tcga.cli=crbind2DataFrame(tcga.cli)
tcga.hall.cox=cox_batch(tcga.hall.ssGSEA[,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
table(tcga.hall.cox$p.value<0.05)
tcga.hall.cox.fit=tcga.hall.cox[tcga.hall.cox$p.value<0.05,]



GSE20685.cli=crbind2DataFrame(GSE20685.cli)
GSE20685.hall.cox=cox_batch(GSE20685.hall.ssGSEA[,GSE20685.cli$Samples],time = GSE20685.cli$OS.time,event = GSE20685.cli$OS)
table(GSE20685.hall.cox$p.value<0.05)
GSE20685.hall.cox.fit=GSE20685.hall.cox[GSE20685.hall.cox$p.value<0.05,]

intersect(rownames(tcga.hall.cox.fit),rownames(GSE20685.hall.cox.fit))
GSE20685.hall.cox.fit[intersect(rownames(tcga.hall.cox.fit),rownames(GSE20685.hall.cox.fit)),]
tcga.hall.cox.fit[intersect(rownames(tcga.hall.cox.fit),rownames(GSE20685.hall.cox.fit)),]


tcga.cli.merge=data.frame(tcga.cli,score=tcga.hall.ssGSEA['HALLMARK_GLYCOLYSIS',tcga.cli$Samples])
tcga.cli.merge=crbind2DataFrame(tcga.cli.merge)
head(tcga.cli.merge)
summary(coxph(formula=Surv(OS.time, OS)~score, data=tcga.cli.merge))
summary(coxph(formula=Surv(DSS.time, DSS)~score, data=tcga.cli.merge))
summary(coxph(formula=Surv(PFS.time, PFS)~score, data=tcga.cli.merge))
library(survival)
library(survminer)
cutoff<-surv_cutpoint(tcga.cli.merge,
                      time="DSS.time",
                      event="DSS",
                      variables='score')
summary(cutoff)
tcga.cli.merge$group=ifelse(tcga.cli.merge$score>cutoff$cutpoint$cutpoint,'High','Low')
tcga.glycolysis.km=ggsurvplot(fit=survfit( Surv(DSS.time/365, DSS) ~ group,
                                       data = tcga.cli.merge),
                          data=tcga.cli.merge,
                          conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                          linetype = c("solid", "dashed","strata")[1],
                          legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                          legend.title = "group",palette = pal_d3()(10)[c(7,9)],
                          legend.labs = c("High","Low"))+ggtitle('TCGA(DSS)')
tcga.glycolysis.km=mg_merge_plot(tcga.glycolysis.km$plot,tcga.glycolysis.km$table,nrow=2,heights = c(3,1),align = 'v')
tcga.glycolysis.km
ggsave('results/01.prognosis/TCGA_glycolysis_km_DSS.pdf',tcga.glycolysis.km,height = 5,width = 4)




GSE20685.cli.merge=data.frame(GSE20685.cli,score=GSE20685.hall.ssGSEA['HALLMARK_GLYCOLYSIS',GSE20685.cli$Samples])
GSE20685.cli.merge=crbind2DataFrame(GSE20685.cli.merge)
head(GSE20685.cli.merge)
summary(coxph(formula=Surv(OS.time, OS)~score, data=GSE20685.cli.merge))
cutoff<-surv_cutpoint(GSE20685.cli.merge,
                      time="OS.time",
                      event="OS",
                      variables='score')
summary(cutoff)
GSE20685.cli.merge$group=ifelse(GSE20685.cli.merge$score>cutoff$cutpoint$cutpoint,'High','Low')
GSE20685.glycolysis.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ group,
                                           data = GSE20685.cli.merge),
                              data=GSE20685.cli.merge,
                              conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                              # title='Stage I+II',
                              linetype = c("solid", "dashed","strata")[1],
                              legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                              legend.title = "group",palette = pal_d3()(10)[c(7,9)],
                              legend.labs = c("High","Low"))+ggtitle('GSE20685(OS)')
GSE20685.glycolysis.km=mg_merge_plot(GSE20685.glycolysis.km$plot,GSE20685.glycolysis.km$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE20685.glycolysis.km
ggsave('results/01.prognosis/GSE20685.glycolysis.km.pdf',GSE20685.glycolysis.km,height = 5,width = 4)

#####1.2 糖酵解与致癌信号通路############
library(progeny)
tcga.pathway.activ=progeny(as.matrix(tcga.exp),scale = T)
tcga.pathway.activ[1:5,1:5]

df=cbind.data.frame(GLYCOLYSIS=tcga.hall.ssGSEA['HALLMARK_GLYCOLYSIS',tcga.cli$Samples],
                    tcga.pathway.activ[tcga.cli$Samples,])
pdf('results/01.prognosis/GLYCOLYSIS_cor_path.pdf',height = 6,width = 6)
corrplot(cor(df),
         tl.col = 'black',
         # cl.align.text = 'l',
         tl.srt=90, diag = F,
         col=colorRampPalette(c('blue', 'white','red'))(100),
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[7],
         order = c("original", "AOE", "FPC", "hclust", "alphabet")[1],
         type="upper")
dev.off()



#######1.3 糖酵解与免疫############
library(ggcorrplot)
library(psych)
tcga.mcp=immu_MCPcounter(tcga.exp)
glycolysis_cor_immu <- corr.test(x =tcga.hall.ssGSEA['HALLMARK_GLYCOLYSIS',tcga.cli$Samples],
                         y = tcga.mcp[tcga.cli$Samples,],
                         method = "spearman",adjust = "BH",ci = F)
glycolysis_cor_immu_res=data.frame(Immune=colnames(tcga.mcp))
glycolysis_cor_immu_res$cor<-as.numeric(glycolysis_cor_immu$r)
glycolysis_cor_immu_res$p.adj<-as.numeric(glycolysis_cor_immu$p.adj)
head(glycolysis_cor_immu_res)
table(glycolysis_cor_immu_res$p.adj<0.05)
glycolysis_cor_immu_res=glycolysis_cor_immu_res[order(glycolysis_cor_immu_res$cor),]
head(glycolysis_cor_immu_res)
library(rcartocolor)
ggplot(data=glycolysis_cor_immu_res,aes(x=cor,y=reorder(Immune,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_viridis_c(option = "plasma") +
  geom_segment(aes(yend=Immune,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Immune cell')+theme_bw()+
  theme(text = element_text(family = 'Times'))
ggsave('results/01.prognosis/GLYCOLYSIS_cor_immu1.pdf',height = 6,width = 6)


tcga.timer=immu_timer(tcga.exp)
glycolysis_cor_immu <- corr.test(x =tcga.hall.ssGSEA['HALLMARK_GLYCOLYSIS',tcga.cli$Samples],
                                 y = tcga.timer[tcga.cli$Samples,],
                                 method = "spearman",adjust = "BH",ci = F)
glycolysis_cor_immu_res=data.frame(Immune=colnames(tcga.timer))
glycolysis_cor_immu_res$cor<-as.numeric(glycolysis_cor_immu$r)
glycolysis_cor_immu_res$p.adj<-as.numeric(glycolysis_cor_immu$p.adj)
head(glycolysis_cor_immu_res)
table(glycolysis_cor_immu_res$p.adj<0.05)
glycolysis_cor_immu_res=glycolysis_cor_immu_res[order(glycolysis_cor_immu_res$cor),]
head(glycolysis_cor_immu_res)
library(rcartocolor)
ggplot(data=glycolysis_cor_immu_res,aes(x=cor,y=reorder(Immune,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_viridis_c(option = "plasma") +
  geom_segment(aes(yend=Immune,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Immune cell')+theme_bw()+
  theme(text = element_text(family = 'Times'))
ggsave('results/01.prognosis/GLYCOLYSIS_cor_immu2.pdf',height = 6,width = 6)


#####02.差异基因#########
dir.create('results/02.DEGs')
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#确定点的颜色
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),#不显示图例标题
      legend.position = leg.pos,
      text = element_text(family = 'Times')
    )+
    ylab(ylab)+#修改y轴名称
    xlab(xlab)+#修改x轴名称
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
  return(p)
}
tcga.limma=mg_limma_DEG(exp = tcga.exp[intersect(mrna_genecode$SYMBOL,rownames(tcga.exp)),tcga.cli.merge$Samples],
                        group = tcga.cli.merge$group, ulab = 'High',dlab = 'Low')
tcga.limma$Summary
tcga.degs=tcga.limma$DEG[tcga.limma$DEG$adj.P.Val<0.05 & abs(tcga.limma$DEG$logFC) > log2(1.5),]
write.csv(tcga.degs,'results/02.DEGs/tcga.degs.csv')
tcga.volcano=my_volcano(dat = tcga.limma,p_cutoff = 0.05,fc_cutoff = log2(1.5))+ggtitle('TCGA')

GSE20685.limma=mg_limma_DEG(exp = GSE20685.exp[intersect(mrna_genecode$SYMBOL,rownames(GSE20685.exp)),GSE20685.cli.merge$Samples],
                            group = GSE20685.cli.merge$group, ulab = 'High',dlab = 'Low')
GSE20685.limma$Summary
GSE20685.degs=GSE20685.limma$DEG[GSE20685.limma$DEG$adj.P.Val<0.05 & abs(GSE20685.limma$DEG$logFC) > log2(1.5),]
dim(GSE20685.degs)
write.csv(GSE20685.degs,'results/02.DEGs/GSE20685.degs.csv')
GSE20685.volcano=my_volcano(dat = GSE20685.limma,p_cutoff = 0.05,fc_cutoff = log2(1.5))+ggtitle('GSE20685')


degs=intersect(rownames(tcga.degs),rownames(GSE20685.degs))
length(degs)
library(eulerr)
v=list(TCGA=rownames(tcga.degs),GSE20685=rownames(GSE20685.degs))
# names(v)=c('CRD related genes','TCGA DEGs','HCCDB18 DEGs')
venn.plot=plot(venn(v),labels = list(col = "gray20", font = 2), 
           edges = list(col="gray60", lex=1),
           fills = list(fill = c("#297CA0", "#E9EA77"), alpha = 0.6),
           quantities = list(cex=.8, col='gray20'))
venn.plot

degs.enrichmenr=mg_clusterProfiler(degs)
p=enrichplot::dotplot(degs.enrichmenr$KEGG)+theme(text = element_text(family = 'Times'))
write.csv(degs.enrichmenr$KEGG@result,'results/02.DEGs/DEGs_enrichment.csv',row.names = F)

fig2=mg_merge_plot(tcga.volcano,GSE20685.volcano,venn.plot,p,ncol=2,nrow=2,labels = LETTERS[1:4])
ggsave('results/02.DEGs/Fig2.pdf',fig2,height = 10,width = 14)

dotplot(degs.enrichmenr$GO_BP)+theme(text = element_text(family = 'Times'))
#######03.风险模型#############
dir.create('results/03.model')
pre.gene=intersect(rownames(tcga.degs),rownames(GSE20685.degs))
length(pre.gene)
#299

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.gene, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

tcga.cli=crbind2DataFrame(tcga.cli)
gene.cox=cox_batch(dat = tcga.exp[pre.gene,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
gene.cox=na.omit(gene.cox)
table(gene.cox$p.value<0.01)
gene.cox=gene.cox[gene.cox$p.value<0.01,]
rownames(gene.cox)= gsub('-', '_',rownames(gene.cox))
gene.cox
write.csv(gene.cox,'results/03.model/gene_cox.csv')

setdiff(rownames(gene.cox),colnames(tcga_model_data))
library(glmnet)
set.seed(2023)
fit1=glmnet(as.matrix(tcga_model_data[,rownames(gene.cox)])
            ,cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix( tcga_model_data[,rownames(gene.cox)])
                  ,cbind(time=tcga_model_data$OS.time,
                         status=tcga_model_data$OS)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)
length(names(sig.coef))

pdf('results/03.model/LASSO.pdf',height = 5,width = 10)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

#######逐步回归
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)
lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
ggforest(cox, data = tcga_model_data, 
         main = "Hazardratio", fontsize =1.0, 
         noDigits = 2)+theme(text = element_text(family = 'Times'))
ggsave('results/03.model/module.coxforest.pdf',height = 5,width = 7)

#############TCGA#####
risktype.col=pal_d3()(10)[c(7,9)]
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
tcga.risktype.cli=crbind2DataFrame(tcga.risktype.cli)
tcga.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
                            variables = 'Riskscore')
tcga.point.cutoff <- as.numeric(summary(tcga.point)[1])
tcga.point.cutoff

tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.point.cutoff,'High','Low')
tcga.roc.OS=ggplotTimeROC(tcga.risktype.cli$OS.time,
                          tcga.risktype.cli$OS,
                          tcga.risktype.cli$Riskscore,mks = c(1,3,5))
tcga.roc.OS
# ggsave('results/04.model/tcga.roc.OS.pdf',height = 5,width = 5)

tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,#ylab='Overall Survival(OS)',
                      legend.position='top',
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.OS


######GSE20685##########
GSE20685_model_data <- cbind(GSE20685.cli[, c("OS.time", "OS")],
                             t(GSE20685.exp[pre.gene, GSE20685.cli$Samples]))
colnames(GSE20685_model_data) <- gsub('-', '_', colnames(GSE20685_model_data))


setdiff(names(lan),colnames(GSE20685_model_data))

risk.GSE20685=as.numeric(lan%*%as.matrix(t(GSE20685_model_data[GSE20685.cli$Samples,names(lan)])))
GSE20685.risktype.cli=data.frame(GSE20685.cli,Riskscore=risk.GSE20685)
GSE20685.point <- surv_cutpoint(GSE20685.risktype.cli, time = "OS.time", event = "OS",
                                variables = 'Riskscore')
GSE20685.point.cutoff <- as.numeric(summary(GSE20685.point)[1])
GSE20685.point.cutoff

GSE20685.risktype.cli$Risktype=ifelse(GSE20685.risktype.cli$Riskscore>GSE20685.point.cutoff,'High','Low')
GSE20685.roc.OS=ggplotTimeROC(GSE20685.risktype.cli$OS.time,
                              GSE20685.risktype.cli$OS,
                              GSE20685.risktype.cli$Riskscore,mks = c(1,3,5))
GSE20685.roc.OS
GSE20685.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                      data = GSE20685.risktype.cli),
                          data=GSE20685.risktype.cli,
                          conf.int = T,pval = T,risk.table = T, 
                          fun = "pct",size = 1,surv.median.line = 'hv',
                          title='GSE20685',legend.title='Risktype',
                          legend.labs = c('High','Low'),
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,#ylab='Overall Survival(OS)',
                          legend.position='top',
                          ggtheme = theme_bw(base_size = 12))
GSE20685.km.OS=mg_merge_plot(GSE20685.km.OS$plot,GSE20685.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE20685.km.OS


roc_km_plot=mg_merge_plot(tcga.km.OS,tcga.roc.OS,GSE20685.km.OS,GSE20685.roc.OS,ncol=2,nrow=2,labels = LETTERS[4:7])
ggsave('results/03.model/roc_km_plot.pdf',roc_km_plot,height = 10,width = 10)

my_mutiboxplot( t(tcga.exp[names(lan),tcga.risktype.cli$Samples]),notch = T,
                group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Gene Expression',fill = 'Risktype',angle = 0,hjust = .5,title = 'TCGA')
ggsave('results/03.model/TCGA_expr.pdf',height = 5,width = 7)
my_mutiboxplot( t(GSE20685.exp[names(lan),GSE20685.risktype.cli$Samples]),
                group = GSE20685.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Gene Expression',fill = 'Risktype',angle = 0,hjust = .5,notch = T,title = 'GSE20685')
ggsave('results/03.model/GSE20685_expr.pdf',height = 5,width = 7)

library(ggbiplot)
tcga.pca <- prcomp(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]), scale=T)
pca_fig1 <- ggbiplot(tcga.pca, scale=1, groups = tcga.risktype.cli$Risktype,
                     ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values =risktype.col) + 
  theme(legend.direction = 'horizontal', legend.position = 'top',text = element_text(family = 'Times')) +
  xlab('PCA1') + ylab('PCA2')+ggtitle('TCGA')
pca_fig1

GSE20685.pca <- prcomp(t(GSE20685.exp[names(lan),GSE20685.risktype.cli$Samples]), scale=T)
pca_fig2 <- ggbiplot(GSE20685.pca, scale=1, groups = GSE20685.risktype.cli$Risktype,
                     ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = risktype.col) + 
  theme(legend.direction = 'horizontal', legend.position = 'top',text = element_text(family = 'Times')) +
  xlab('PCA1') + ylab('PCA2')+ggtitle('GSE20685')
pca_fig2


mg_merge_plot(pca_fig1,pca_fig2,labels = c('I','J'),align = 'h',common.legend = T)
ggsave('results/03.model/pca.pdf',height = 5,width = 10)



#######04.临床特征########
dir.create('results/04.clinical')
p1=my_violin(dat = tcga.risktype.cli$Riskscore,group = tcga.risktype.cli$T.stage,group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'T stage')
p2=my_violin(dat = tcga.risktype.cli$Riskscore,group = tcga.risktype.cli$N.stage,group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'N stage')
p3=my_violin(dat = tcga.risktype.cli$Riskscore,group = tcga.risktype.cli$M.stage,group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'M stage')
p4=my_violin(dat = tcga.risktype.cli$Riskscore,group = tcga.risktype.cli$Stage,group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'Stage')



p5=my_violin(dat = GSE20685.risktype.cli$Riskscore[GSE20685.risktype.cli$T.stage!='TNA'],
          group = GSE20685.risktype.cli$T.stage[GSE20685.risktype.cli$T.stage!='TNA'],group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'T stage')
p6=my_violin(dat = GSE20685.risktype.cli$Riskscore,
          group = GSE20685.risktype.cli$N.stage,group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'N stage')
p7=my_violin(dat = GSE20685.risktype.cli$Riskscore,group = GSE20685.risktype.cli$M.stage,group_cols = pal_jco()(10),
          ylab = 'Riskscore',legend.position = 'none',xlab = 'M stage')

clinical.plot=mg_merge_plot(mg_merge_plot(p1,p2,p3,p4,ncol=4,widths = c(1,1,0.6,1)),
              mg_merge_plot(p5,p6,p7,ncol=3,widths = c(1,1,0.6)),
              nrow=2,labels = c('A','B'))
ggsave('results/04.clinical/clinical.plot.pdf',clinical.plot,height = 9,width = 15)


tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'


table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)

#######单因素#####
#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat


#AJCC_stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <- c("Age",
                        "T.stage",
                        "N.stage",
                        "M.stage",
                        "Stage",
                        "RiskScore")
data.sig$Features=rownames(data.sig) 
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
write.xlsx(data.sig,'results/04.clinical/Univariate_cox.xlsx',overwrite = T)
pdf('results/04.clinical/Univariate.pdf',height = 5,width = 7,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap = 8,lineheight = 10,graph.pos = 3)
dev.off()


#多因素
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+T.stage +  N.stage +M.stage+Stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Features=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Features,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age",
                         "T.stage",
                         "N.stage",
                         "M.stage",
                         "Stage",
                         "RiskScore")
data.muti$Features <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
write.xlsx(data.muti,'results/04.clinical/Multivariate_cox.xlsx',overwrite = T)
pdf('results/04.clinical/Multivariate.pdf',height = 5,width = 7,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap = 8,lineheight = 10,graph.pos = 3)
dev.off()

######列线图
pdf('results/04.clinical/nomogram.pdf', width = 12, height = 10)
mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
                     quick=T,mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#对观测2的六个指标在列线图上进行计分展示
  #,observation=pbc[2,] #也可以不展示
  #预测3年和5年的死亡风险，此处单位是day
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox回归中需要TRUE
  #              ,showP = T #是否展示统计学差异
  #              ,droplines = F#观测2示例计分是否画线
  #,colors = mg_colors[1:3] #用前面自己定义的颜色
  #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #展示观测的可信区间
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                Stage=tcga_cox_datas$Stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))

########05.TME##########
dir.create('results/05.TME')
######TCGA免疫浸润
tcga.estimate=read.delim('results/05.TME/TCGA_ESTIMATE_score.txt')
head(tcga.estimate)
tcga.immu.p1=my_mutiboxplot( tcga.estimate[tcga.risktype.cli$Samples,1:3],notch = T,
                   group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                   ylab = 'Score',fill = 'Risktype',angle = 30)

tcga.mcp=immu_MCPcounter(tcga.exp)
tcga.immu.p2=my_mutiboxplot( tcga.mcp[tcga.risktype.cli$Samples,],notch = T,
                group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Score',fill = 'Risktype',angle = 30)

tcga.cibersort=read.delim('results/05.TME/tcga_cibersort.txt',row.names = 1,check.names = F)
head(tcga.cibersort)
tcga.immu.p3=my_mutiboxplot(tcga.cibersort[tcga.risktype.cli$Samples,1:22],notch = T,
                             group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                             ylab = 'Score',fill = 'Risktype',angle = 30)

#########GSE20685免疫浸润
GSE20685.estimate=read.delim('results/05.TME/GSE20685_ESTIMATE_score.txt')
head(GSE20685.estimate)
GSE20685.immu.p1=my_mutiboxplot( GSE20685.estimate[GSE20685.risktype.cli$Samples,1:3],notch = T,
                group = GSE20685.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Score',fill = 'Risktype',angle = 30)

GSE20685.mcp=immu_MCPcounter(GSE20685.exp)
GSE20685.immu.p2=my_mutiboxplot(GSE20685.mcp[GSE20685.risktype.cli$Samples,],notch = T,
                group = GSE20685.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Score',fill = 'Risktype',angle = 30)


GSE20685.cibersort=read.delim('results/GSE20685_cibersort.txt',row.names = 1,check.names = F)
head(GSE20685.cibersort)
GSE20685.immu.p3=my_mutiboxplot(GSE20685.cibersort[GSE20685.risktype.cli$Samples,1:22],notch = T,
                group = GSE20685.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Score',fill = 'Risktype',angle = 30)



immune.plot=mg_merge_plot(mg_merge_plot(tcga.immu.p1,tcga.immu.p2,tcga.immu.p3,ncol=3,widths = c(1.5,3,5),
                                        align = 'h',common.legend = T,labels = LETTERS[1:3]),
                          mg_merge_plot(GSE20685.immu.p1,GSE20685.immu.p2,GSE20685.immu.p3,ncol=3,widths = c(1.5,3,5),
                                        align = 'h',common.legend = T,labels = LETTERS[4:6]),
                          nrow=2,common.legend = T)
ggsave('results/05.TME/immune.plot.pdf',immune.plot,height = 10,width = 20)


#####06.免疫治疗########
dir.create('results/06.immune_treatment')
####免疫表观评分
tcia.res=read.delim('results/06.immune_treatment/TCIA-ClinicalData.tsv')
tcia.res=data.frame(Samples=paste0(tcia.res$barcode,'-01'),
                    tcia.res[,26:29])
rownames(tcia.res)=tcia.res$Samples
tcia.res=tcia.res[,-1]
colnames(tcia.res)=c('IPS','IPS-CTLA4 blocker','IPS-PD1/PDL1/PDL2 blocker','IPS-CTLA4-and \nPD1/PDL1/PDL2 blocker')
head(tcia.res)
table(tcia.res$IPS)
table(tcia.res$`IPS-CTLA4 blocker`)
table(tcia.res$`IPS-PD1/PDL1/PDL2 blocker`)




my_mutiboxplot(tcia.res[tcga.risktype.cli$Samples,],
                group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Score',fill = 'Risktype',angle = 0,hjust = .5)
ggsave('results/06.immune_treatment/IPS.pdf',height = 5,width = 7)

####抗原呈递
immunomodulators=readxl::read_xlsx('Immunomodulators_PMID29628290.xlsx')
immunomodulators=immunomodulators[,c(2,6)]
immunomodulators.gene=immunomodulators[immunomodulators$`HGNC Symbol` %in% rownames(tcga.exp),]
immunomodulators.gene=crbind2DataFrame(immunomodulators.gene)
immunomodulators.gene=immunomodulators.gene[order(immunomodulators.gene$`Super Category`),]

pdf('results/06.immune_treatment/Antigen_presentation.pdf',height = 5,width = 6)
Heatmap(as.matrix(t(scale(t(tcga.exp[immunomodulators.gene$`HGNC Symbol`[immunomodulators.gene$`Super Category`=='Antigen presentation'],
                                     tcga.risktype.cli$Samples[order(tcga.risktype.cli$Riskscore)]])))),
        name = "Expr",
        column_split = tcga.risktype.cli$Risktype[order(tcga.risktype.cli$Riskscore)],
        column_title_gp = gpar(fill =risktype.col),
        cluster_rows = F, cluster_columns = F,
        cluster_row_slices = F, cluster_column_slices=F,
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF')))
dev.off()


##########免疫检查点抑制剂
icgs.gene=read.xlsx('79ICGs_PMID32814346.xlsx',check.names = F)
head(icgs.gene)
icgs.gene=icgs.gene[,c(1,3)]
table(icgs.gene$Role.with.Immunity)
icgs.gene=icgs.gene[icgs.gene$Symbol%in%rownames(tcga.exp),]
icgs.gene=icgs.gene[icgs.gene$Role.with.Immunity!='TwoSide',]
my_mutiboxplot(t(tcga.exp[icgs.gene$Symbol[icgs.gene$Role.with.Immunity=='Inhibit'],tcga.risktype.cli$Samples]),
               group = tcga.risktype.cli$Risktype,legend.pos = 'none',group_cols = risktype.col,
               ylab = 'Gene Expression',notch = T)
ggsave('results/06.immune_treatment/ICGs.pdf',height = 5,width = 14)

###########药物敏感性
load('results/06.immune_treatment/tcga_durg_ic50.RData')
tcga_durg_ic50_res[1:5,1:5]
colnames(tcga_durg_ic50_res)[1]='Cisplatin'

IC50_cor_RS <- corr.test(x = tcga.risktype.cli$Riskscore,
                         y = tcga_durg_ic50_res[tcga.risktype.cli$Samples,],
                         method = "spearman",adjust = "BH",ci = F)
IC50_cor_RS_res=data.frame(Drugs=colnames(tcga_durg_ic50_res))
IC50_cor_RS_res$cor<-as.numeric(IC50_cor_RS$r)
IC50_cor_RS_res$p.adj<-as.numeric(IC50_cor_RS$p.adj)
head(IC50_cor_RS_res)
table(IC50_cor_RS_res$p.adj<0.05 , abs(IC50_cor_RS_res$cor) > 0.3)
IC50_cor_RS_res=IC50_cor_RS_res[IC50_cor_RS_res$p.adj<0.05 & abs(IC50_cor_RS_res$cor) > 0.3,]
IC50_cor_RS_res=IC50_cor_RS_res[order(IC50_cor_RS_res$cor),]
head(IC50_cor_RS_res)
library(rcartocolor)
ggplot(data=IC50_cor_RS_res,aes(x=cor,y=reorder(Drugs,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_carto_c(palette = "Earth")+
  geom_segment(aes(yend=Drugs,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Drugs')+theme_bw()+
  theme(text = element_text(family = 'Times'))

ggsave('results/06.immune_treatment/IC50_cor_RS.pdf',height = 5,width = 5)


my_mutiboxplot(tcga_durg_ic50_res[tcga.risktype.cli$Samples,IC50_cor_RS_res$Drugs],
               group = tcga.risktype.cli$Risktype,legend.pos = 'none',group_cols = risktype.col,
               ylab = 'IC50',notch = T)
ggsave('results/06.immune_treatment/drug_sensitive.pdf',height = 5,width = 9)

save.image(file='project.RData')


