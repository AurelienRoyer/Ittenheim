### Package library needed
library(devtools)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(vegan)
### 95%  confidence intervals binomial
bt <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=c("two.sided"), conf.level = 0.95)$conf.int[1:2]}

### Load data
Data.tapho<-read.csv("supplementary data.csv", header = T,sep=";",row.names="Code_uniq")
prey<-as.factor(Data.tapho$Prey)
MNI<-Data.tapho$MNI
cat<-Data.tapho$Type.of.accumulation

# Cleaning data used
bones.selected<-c("MANDIBLE","FEMUR","TIBIA","PELVIS","HUMERUS","RADIUS","ULNA","SCAPULA")
Data.tapho.sq<-Data.tapho[,c(bones.selected)]

### Calculation of ratios
CRA.POSTCRA<-Data.tapho.sq$MAND*2/(Data.tapho.sq$MAND*2+Data.tapho.sq$HUM+Data.tapho.sq$FEM)
t.CRA.POSTCRA<-data.frame(t(mapply(bt,Data.tapho.sq$MAND*2,Data.tapho.sq$MAND*2+Data.tapho.sq$HUM+Data.tapho.sq$FEM)),CRA.POSTCRA)

AN.PO<-(Data.tapho.sq$HUM+Data.tapho.sq$ULNA) / (Data.tapho.sq$HUM+Data.tapho.sq$ULNA+Data.tapho.sq$FEM+Data.tapho.sq$TIB) 
t.AN.PO<-data.frame(t(mapply(bt,Data.tapho.sq$HUM+Data.tapho.sq$ULNA,Data.tapho.sq$HUM+Data.tapho.sq$ULNA+Data.tapho.sq$FEM+Data.tapho.sq$TIB)),AN.PO)

Z.E<-(Data.tapho.sq$TIB+Data.tapho.sq$ULNA) / (Data.tapho.sq$TIB+Data.tapho.sq$ULNA+Data.tapho.sq$FEM+Data.tapho.sq$HUM)
t.Z.E<-data.frame(t(mapply(bt,Data.tapho.sq$TIB+Data.tapho.sq$ULNA,Data.tapho.sq$TIB+Data.tapho.sq$ULNA+Data.tapho.sq$FEM+Data.tapho.sq$HUM)),Z.E)
ALL.ratio<-data.frame(t.CRA.POSTCRA*100,t.AN.PO*100,t.Z.E*100)

### missMDA
ncomp<-estim_ncpPCA(Data.tapho.sq)
print(ncomp$ncp)
res.imput<-imputeCA(Data.tapho.sq, ncp=ncomp$ncp)

### Canonical analyses
res2<-CA(res.imput, row.sup = c(1:2))

#number of axes to keep
res2$eig
fviz_screeplot(res2, addlabels = TRUE, ylim = c(0, 50))+
  geom_hline (yintercept = 1/ncol(res.imput-1)*100, linetype = 2, color = "red")

#contribution of dim 1, 2 and 3
res2$col$contrib
fviz_contrib(res2, choice = "col", axes = 1)
fviz_contrib(res2, choice = "col", axes = 2)
fviz_contrib(res2, choice = "col", axes = 3)
fviz_contrib (res2, choice = "col", axes = 1:3)


### Groups for shape, color and size
shape.test <- as.numeric(as.character(fct_recode(prey,"23" = "Talpa","21" = "Rodent","22" = "Leporid")))
color.test<- as.character(fct_recode(cat,"#CC0099" = "Den","#000000" = "Ittenheim","#56B4E9" = "Mix",
                                     "#CC3300" = "Natural death","#009966" = "Non-ingested","#999999" = "Pellets","#FF9900" = "Scat"))
MNI2<-case_when(MNI<=15 ~2,MNI>100~10,MNI>11~5)

### Figures CA analysis

p<-fviz_ca_biplot(res2, axes = c(1,2), col.row.sup = color.test[1:2],  col.row=as.factor(cat[3:length(cat)]), repel=TRUE,label="all",mean.point=0)
p+geom_point(fill=color.test[3:length(color.test)],size=MNI2[3:length(MNI2)],shape=shape.test[3:length(shape.test)])+
  scale_color_manual(values=c("#FF3399","#56B4E9","red","#009966","#999999","#FF9900"))

p<-fviz_ca_biplot(res2, axes = c(1,3), col.row.sup = color.test[1:2], col.row=as.factor(cat[3:length(shape.test)]), repel=TRUE,label="all",mean.point=0)
p+geom_point(fill=color.test[3:length(color.test)],size=MNI2[3:length(MNI2)],shape=shape.test[3:length(shape.test)])+
  scale_color_manual(values=c("#FF3399","#56B4E9","red","#009966","#999999","#FF9900"))


### Figures Plot

ggplot(data=ALL.ratio, aes(x=ALL.ratio$CRA.POSTCRA,y=ALL.ratio$ AN.PO),label=cat)+
  geom_errorbar(aes(ymin=X1.1, ymax=X2.1),col="grey",alpha=0.25)+
  geom_errorbarh(aes(xmin=X1, xmax=X2),col="grey",alpha=0.25)+
  geom_point(fill=color.test,size=MNI2,shape=shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames(Data.tapho)),data=Data.tapho,size=3)+
  xlab("CRA/POSTCRA%")+ylab("AN/PO%")+
  theme(legend.position="bottom")

ggplot(data=ALL.ratio, aes(x=ALL.ratio$CRA.POSTCRA,y=ALL.ratio$Z.E),label=cat)+
  geom_errorbar(aes(ymin=X1.2, ymax=X2.2),col="grey",alpha=0.25)+
  geom_errorbarh(aes(xmin=X1, xmax=X2),col="grey",alpha=0.25)+
  geom_point(fill=color.test,size=MNI2,shape=shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames(Data.tapho)),data=Data.tapho,size=3)+
  xlab("CRA/POSTCRA%")+ylab("Z/E%")+
  theme(legend.position="bottom")

#without error bars
ggplot(data=Data.tapho, aes(x=cra.postcra*100,y=AN.PO*100),label=cat)+
  geom_point(fill=color.test,size=MNI2,shape=shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames(Data.tapho)),data=Data.tapho,size=3)+
  xlab("CRA/POSTCRA%")+ylab("AN/PO%")+
  theme(legend.position="bottom")

ggplot(data=Data.tapho, aes(x=cra.postcra*100,y=Z.E*100),label=cat)+
  geom_point(fill=color.test,size=MNI2,shape=shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames(Data.tapho)),data=Data.tapho,size=3)+
  xlab("CRA/POSTCRA%")+ylab("Z/E%")+
  theme(legend.position="bottom")


### Correlation between CA dimensions and skeletal ratios
cor.test(as.numeric(as.character(Data.tapho$cra.postcra[3:nrow(Data.tapho)])),res2$row$coord[,1])
cor.test(as.numeric(as.character(Data.tapho$AN.PO[3:nrow(Data.tapho)])),res2$row$coord[,2])
cor.test(as.numeric(as.character(Data.tapho$Z.E[3:nrow(Data.tapho)])),res2$row$coord[,3])


###Figure plots digestion
#digested femur vs incisors
#data preparation
data.cat<-data.frame(x1=c(0,20,30,50,80),x2=c(20,40,50,80,100),x3=c(0,20,50,60,80), x4=c(15,30,70,80,100),r=c("cat1","cat2","cat3","cat4","cat5"))
Data.digIfem<- select(Data.tapho,(c("FEMUR","INCISORS","X.Femur.digested", "X.Incisors.digested","X..tooth.marks"))) %>% 
                                       data.frame(rownames(Data.tapho),cat,color.test,MNI2,shape.test,rowSums(Data.tapho.sq[2:7]))%>%
                                    filter(!is.na(Data.tapho$X.Incisors.digested) & !is.na(Data.tapho$X.Femur.digested))
t.dig<-data.frame(Data.digIfem,
           t(mapply(bt,round(Data.digIfem$X.Incisors.digested*Data.digIfem$INCISORS/100,0),Data.digIfem$INCISORS,Data.digIfem)*100),
           t(mapply(bt,round(Data.digIfem$X.Femur.digested*Data.digIfem$FEMUR/100,0),Data.digIfem$FEMUR,Data.digIfem)*100))

#figure
ggplot(data=t.dig, aes(x=X.Femur.digested, y=X.Incisors.digested),label=cat)+
  geom_errorbar(aes(ymin=X1, ymax=X2),col="grey",alpha=0.25)+
  geom_errorbarh(aes(xmin=X1.1, xmax=X2.1),col="grey",alpha=0.25)+
  geom_point(fill=t.dig$color.test,size=t.dig$MNI2,shape=t.dig$shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames.Data.tapho.),data=t.dig,size=3)+
  geom_rect(data=data.cat,aes(xmin=x1,xmax=x2,ymin=x3,ymax=x4),fill= "red", alpha=0.2, inherit.aes = FALSE) +
  geom_text(data=data.cat, aes(x=x1, y=x4+1, label=r), col="blue")+
  xlab("%Femur digestion")+ylab("%Incisor digestion")+
  theme(legend.position="bottom")


#without error bars
ggplot(data=Data.tapho, aes(x=X.Femur.digested, y=X.Incisors.digested),label=cat)+
  geom_point(fill=color.test,size=MNI2,shape=shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=row.names(Data.tapho)),data=Data.tapho,size=3)+
  geom_rect(data=data.cat,aes(xmin=x1,xmax=x2,ymin=x3,ymax=x4),fill= "red", alpha=0.2, inherit.aes = FALSE) +
  geom_text(data=data.cat, aes(x=x1, y=x4+1, label=r), col="blue",inherit.aes = FALSE)+
  xlab("%Femur digestion")+ylab("%Incisor digestion")+
  theme(legend.position="bottom")


# %digested incisor vs toothmarks%
#Data preparation
Data.digIfem2<-filter(t.dig, cat != "Pellets") %>%na.omit(t.dig)
t.tooth<-data.frame(t(mapply(bt,round(Data.digIfem2$X..tooth.marks*Data.digIfem2$rowSums.Data.tapho.sq.2.7../100,0),Data.digIfem2$rowSums.Data.tapho.sq.2.7..,Data.digIfem2))*100,Data.digIfem2)

#Figure
ggplot(data=t.tooth, aes(x=X..tooth.marks, y=X.Incisors.digested),label=cat)+
  geom_errorbar(aes(ymin=X1.1, ymax=X2.1),col="grey",alpha=0.25)+
  geom_errorbarh(aes(xmin=X1.2, xmax=X2.2),col="grey",alpha=0.25)+
  geom_point(fill=t.tooth$color.test,size=t.tooth$MNI2,shape=t.tooth$shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames.Data.tapho.),data=t.tooth,size=3)+
  xlab("%tooth marks")+ylab("%Incisor digestion")+
  theme(legend.position="bottom")

#without error bars
ggplot(data=t.tooth, aes(x=X..tooth.marks, y=X.Incisors.digested),label=cat)+
  geom_point(fill=t.tooth$color.test,size=t.tooth$MNI2,shape=t.tooth$shape.test)+
  scale_color_manual(values=c("#CC0099","#000000","#56B4E9","#CC3300","#009966","#999999","#FF9900"))+
  geom_text_repel(aes(label=rownames.Data.tapho.),data=t.tooth,size=3)+
  xlab("%tooth marks")+ylab("%Incisor digestion")+
  theme(legend.position="bottom")

### NP-MANOVA
#Data preparation
res.imput.idx<-data.frame(res.imput,prey,cat)
res.imput.idx<-  res.imput.idx[Data.tapho$Type.of.accumulation == "Scat" | Data.tapho$Type.of.accumulation == "Natural death"  
                          | Data.tapho$Type.of.accumulation == "Non-ingested" 
                          |Data.tapho$Type.of.accumulation == "Den"
                          &Data.tapho$Prey != "Talpa", ]

#check the dispersion
dis <- vegdist(res.imput.idx[1:(ncol(res.imput.idx)-2)])
mod<-betadisper(dis, res.imput.idx$cat, type = "median", bias.adjust = FALSE,
                sqrt.dist = FALSE, add = FALSE)
anova(mod)
pmod<-permutest(mod, pairwise = TRUE, permutations = 99)
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

mod<-betadisper(dis, res.imput.idx$prey, type = "median", bias.adjust = FALSE,
                sqrt.dist = FALSE, add = FALSE)
anova(mod)
#Perform a one-way NPMANOVA
adonis(res.imput.idx[1:(ncol(res.imput.idx)-2)] ~ res.imput.idx$cat*res.imput.idx$prey,method = "bray",
        per.mutations = 9999)

#post-hoc analysis  
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwise.adonis)
pairwise.adonis(res.imput.idx[1:(ncol(res.imput.idx)-2)], res.imput.idx$cat, sim.method = "bray",
                p.adjust.m = "bonferroni")



#### Code figures and references

print((paste(rownames(Data.tapho), Data.tapho$Reference, sep="-")))

      