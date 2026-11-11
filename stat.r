library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrastr)
library(ggrepel)
library(patchwork)
stat=fread("stat.csv",data.table=F)

df_long <- stat[which(stat$group=="GENE"),c(1,5:7)] %>%
  pivot_longer(cols = c(ncommon, nrawu, nstu), 
               names_to = "category", 
               values_to = "count")

df_percent <- df_long %>%
  group_by(trait) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

df_percent$category=factor(df_percent$category,levels=c("nrawu","nstu","ncommon"))
p=ggplot(df_percent, aes(x = trait, y = percentage, fill = category)) +
  geom_col(position = "fill") + scale_fill_brewer(palette = "Set2") + 
  scale_y_continuous(labels = scales::percent_format())  +
  theme_classic()+ylab("% of significant genes")+theme(
    axis.ticks=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8,color="black"))

pdf("Figure 3A.pdf",height=4,width=4)
p
dev.off()
burden=fread("income.burden",data.table=F)
burden=separate(burden,ID,into=c("gene","type"),sep="\\(")
burden$P=10^(-as.numeric(burden$LOG10P))
burdenst=fread("income.st.burden",data.table=F)
burdenst=separate(burdenst,ID,into=c("gene","type"),sep="\\(")
burdenst$P=10^(-as.numeric(burdenst$LOG10P))

burdensig=burden[which(burden$TEST=="ADD-BURDEN-ACAT"),]
burdensig$padj=p.adjust(burdensig$P,method="fdr")
burdensig[which(burdensig$padj<0.05),]
draw=burdensig[which(burdensig$P<0.05/17705),]

burdenstsig=burdenst[which(burdenst$TEST=="ADD-BURDEN-ACAT"),]
burdenstsig$padj=p.adjust(burdenstsig$P,method="fdr")
burdenstsig[which(burdenstsig$P<0.05/nrow(burdensig)),]


setdiff(burdenstsig[which(burdenstsig$P<0.05/nrow(burdensig)),"gene"],
burdensig[which(burdensig$P<0.05/nrow(burdensig)),"gene"])

gwas_data=burden[which(burden$TEST=="ADD-BURDEN-ACAT"),c("CHROM","GENPOS","P","gene")]
colnames(gwas_data)=c("CHR","BP","P","gene")
gwas_data=gwas_data[-which(gwas_data$gene %in% c("CACNA1E")),]
gwas_data$CHR=as.numeric(gwas_data$CHR)
gwas_data$BP=as.numeric(gwas_data$BP)
gwas_data=gwas_data[!is.na(gwas_data$CHR),]
gwas_data=gwas_data[order(gwas_data$CHR,gwas_data$BP),]
gwas_data$pos=1:nrow(gwas_data)
gwas_data$color_group <- as.factor(gwas_data$CHR %% 2)
gwas_data=gwas_data[!is.na(gwas_data$CHR),]
chr_labels <- gwas_data %>%
  group_by(CHR) %>%
  summarize(chr_mid = mean(pos))
chr_labels =data.frame(chr_labels)


gene=gwas_data[which(gwas_data$P<0.05/17704),]

mraw=ggplot(gwas_data, aes(x = pos, y = -log10(P), color = as.factor(CHR %% 2))) +
  rasterise(geom_point(alpha = 0.8, size = 1),dpi= 720)+
  scale_color_manual(values =c("#47a1a2", "#272973")) +
  scale_x_continuous(breaks = chr_labels$chr_mid[1:22], labels = chr_labels$CHR[1:22],expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0.1))+theme_classic()+
  geom_hline(yintercept=-log10(0.05/17704),color="red")+labs(x ="Chromosome",y = "-log10(P-value)")+
   geom_text_repel(
    data = gene, 
    aes(label = gene),
    size = 4,              # Font size
    max.overlaps = 100,     # Limit the number of overlaps
    box.padding = 0.1,     # Padding around text
    point.padding = 0.1,   # Space around points
    segment.color = "grey" # Line color connecting text to points
  ) + theme(legend.position = "none",axis.line.x=element_blank(),
            axis.text=element_text(color = "black", size = 15),
            ,axis.title=element_text(color = "black", size = 15))
  
gwas_data=burdenst[which(burdenst$TEST=="ADD-BURDEN-ACAT"),c("CHROM","GENPOS","P","gene")]
colnames(gwas_data)=c("CHR","BP","P","gene")
gwas_data$CHR=as.numeric(gwas_data$CHR)
gwas_data$BP=as.numeric(gwas_data$BP)
gwas_data=gwas_data[!is.na(gwas_data$CHR),]
gwas_data=gwas_data[order(gwas_data$CHR,gwas_data$BP),]
gwas_data$pos=1:nrow(gwas_data)
gwas_data$color_group <- as.factor(gwas_data$CHR %% 2)
gwas_data=gwas_data[!is.na(gwas_data$CHR),]
chr_labels <- gwas_data %>%
  group_by(CHR) %>%
  summarize(chr_mid = mean(pos))
chr_labels =data.frame(chr_labels)


gene=gwas_data[which(gwas_data$P<0.05/17704),]

mst=ggplot(gwas_data, aes(x = pos, y = -log10(P), color = as.factor(CHR %% 2))) +
  rasterise(geom_point(alpha = 0.8, size = 1),dpi= 720)+
  scale_color_manual(values =c("#47a1a2", "#272973")) +
  scale_x_continuous(breaks = chr_labels$chr_mid[1:22], labels = chr_labels$CHR[1:22],expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0.1))+theme_classic()+
  geom_hline(yintercept=-log10(0.05/17704),color="red")+labs(x ="Chromosome",y = "-log10(P-value)")+
   geom_text_repel(
    data = gene, 
    aes(label = gene),
    size = 4,              # Font size
    max.overlaps = 100,     # Limit the number of overlaps
    box.padding = 0.1,     # Padding around text
    point.padding = 0.1,   # Space around points
    segment.color = "grey" # Line color connecting text to points
  ) + theme(legend.position = "none",axis.line.x=element_blank(),
            axis.text=element_text(color = "black", size = 15),
            ,axis.title=element_text(color = "black", size = 15))


pdf("Figure 3down.pdf",width=15,height=5)
mraw|mst
dev.off()

library(RColorBrewer)
library(scales)
library(ggpointdensity)
d=fread("annoburden",data.table=F)
b=fread("burden.txt",data.table=F)
d$all=b[match(d[,1],b[,1]),2]
british=fread("British.circadian.pheno",data.table=F)
british=british[,1]
d=d[which(d[,1] %in% british),]


ggplot(d, aes(x = all, y = exon)) +
  rasterise(geom_pointdensity(), dpi = 300)+
  scale_color_gradientn(
    colors = brewer.pal(11, "Spectral")[11:1],  # 使用 Spectral 调色板
    name = "Density")+theme_classic()
ggsave("figure 1scatter.pdf",height=4,width=4)

prop.test(x = sum(d$exon), n = sum(d$all), p = 0.044, correct = TRUE)
