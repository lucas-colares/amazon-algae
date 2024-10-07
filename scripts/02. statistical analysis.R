# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############# Statistical analysis for Palheta et al. (2023) ##############
# # # # # # # # # # By: Lucas Colares & Leandra Palheta # # # # # # # # # #
# # # # # # # # # # # # # # # # 20-08-2023 # # # # # # # ## # # # # # # # #

source("scripts/00. setup.R")

read.csv("datasets/comm.csv",h=T,sep=";",row.names = 1)->comm
read.csv("datasets/env.csv",h=T,sep=";",row.names = 1)->env
read.csv("datasets/HII.csv",h=T,sep=";",row.names = 1)->HII

###### Environment ----
#### Linear models for environmental variables
quimicas=env[,1:9]
fisicas=HII[,1:(ncol(HII)-2)]
landscape=env[,11:15]
# Chemical block
write.csv(cor(scale(data.frame(HII=HII$IIH_total,quimicas))),"results/table s1.csv")
reshape2::melt(cor(scale(data.frame(HII=HII$IIH_total,quimicas,landscape))))->cors
cors[(abs(cors$value)>0.6)&(abs(cors$value)!=1),]

remv<-c("Alk","pH","regeneracao_floresta")
quimicas[,!grepl(paste0(remv,collapse = "|"),colnames(quimicas))]->quimicas
landscape[,!grepl(paste0(remv,collapse = "|"),colnames(landscape))]->landscape

VarNonCor<-data.frame(quimicas[,c("BOD","Phos","Nitg","Amm","Alum","Iron")])
LMData1<-data.frame(scale(VarNonCor),IIH_total=HII$IIH_total)
LM1=lm(IIH_total~.,data=LMData1)
LM1=step(LM1,direction = "both")
summary(LM1)

melt(cor(HII[,-ncol(HII)]))->CorMat
CorMat$p<-NA
for(x in 1:nrow(CorMat)){
  if(CorMat[x,1]==CorMat[x,2]){
    CorMat[x,3]=NA
    next  
  }
  CorTest<-cor.test(HII[,CorMat[x,1]],HII[,CorMat[x,2]])
  CorTest$p.value->CorMat$p[x]
}

CorMat$sig<-NA
CorMat$sig[CorMat$p<=0.001]<-"***"
CorMat$sig[CorMat$p>0.001&CorMat$p<=0.01]<-"**"
CorMat$sig[CorMat$p>0.01&CorMat$p<=0.05]<-"*"
CorMat$sig[CorMat$p>0.05&CorMat$p<=0.1]<-"."

CorMat$value=round(CorMat$value,2)
CorMat[is.na(CorMat$sig),"sig"]<-""
CorMat$value_text=paste0(CorMat$value,CorMat$sig)
CorMat[CorMat$value_text=="NA","value_text"]=""

tiff("raw_figures/Figure 2.tiff", units="in", width=9.5, height=8, res=600)
ggplot(CorMat, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile()+  
  geom_text(aes(label = value_text), color = "black", size = 4)+
  coord_fixed()+
  scale_x_discrete(name="Physical variables",labels=c("Land","Widt","Rip","Veg","Ret","Sed","Bnk","BnkUnd","Bot","Flw","AquVeg","Det","HII"))+
  scale_y_discrete(name="Physical variables",labels=c("Land","Widt","Rip","Veg","Ret","Sed","Bnk","BnkUnd","Bot","Flw","AquVeg","Det","HII"))+
  scale_fill_gradient2(name="Correlation",low = viridis(3)[1],
                       mid = "white",
                       high = viridis(3)[2])+
  theme(legend.text.align = 0,legend.title.align = 0,legend.key.size = unit(0.5,"cm"),legend.position ="right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(family = "sans",size=12, face="bold"), legend.text = element_text(family = "sans",size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))
dev.off()

prcomp(fisicas)->PCA.local
rownames(PCA.local$rotation)=c("LandInt","WidtInt","ComplRip","VegInt","RetInt","SedInt","BnkInt","BnkUnd","BotInt","FlwComp","AquVegInt","DetInt")
write.csv(PCA.local$rotation[,1:2],"results/table s6.csv")

summary(PCA.local)
eigenvals(PCA.local)[1:2] #autovetor. Sempre tem que ser maior que o Broken-Stick
bstick(PCA.local)[1:2] #broken-stick
PCA.local$rotation[,1:2] #loadings. A import?nia de cada vari?vel para os eixos.

PCAloadings1 <- data.frame(Variables = rownames(PCA.local$rotation), PCA.local$rotation)
SelLoad1<-PCAloadings1
SelLoad1[,2][abs(SelLoad1[,2])>0.3]
SelLoad1[,1][abs(SelLoad1[,2])>0.3]
SelLoad1[,3][abs(SelLoad1[,3])>0.3]
SelLoad1[,1][abs(SelLoad1[,3])>0.3]


unique(c(SelLoad1[,1][abs(SelLoad1[,2])>0.3], SelLoad1[,1][abs(SelLoad1[,3])>0.3]))->SelLoads
SelLoad1[grepl(paste0(SelLoads,collapse = "|"),SelLoad1$Variables),]->SelLoad1

pca1<-ggplot(data.frame(PCA.local$x), aes(x = PC1, y = PC2)) +
  geom_hline(yintercept=0, linetype="dashed", alpha=0.5)+
  geom_vline(xintercept=0, linetype="dashed", alpha=0.5)+
  geom_point(aes(x = PC1, y = PC2, color=HII$IIH_total), size=3)+
  geom_segment(data = SelLoad1, aes(x = 0, y = 0, xend = (PC1*1.8),
                                    yend = (PC2*1.8)), arrow = arrow(length = unit(1/4, "picas")),color = "black", size=0.6,alpha=1) +
  geom_text_repel(data=SelLoad1,aes(x = (SelLoad1$PC1*1.8), y = (SelLoad1$PC2*1.8)),
                  label = SelLoad1$Variables, color="black", fontface="italic", size=3,alpha=1)+
  scale_x_continuous(name="PC1 (62.66%)",limits=c(-1.25,1.55))+
  scale_y_continuous(name="PC2 (21.09%)")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  theme(legend.text.align = 0,legend.title.align = 0,legend.key.size = unit(0.5,"cm"),legend.position ="right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(family = "sans",size=12, face="bold"), legend.text = element_text(family = "sans",size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); pca1

prcomp(scale(quimicas))->PCA.quimica
summary(PCA.quimica)
write.csv(PCA.quimica$rotation[,1:2],"results/table s7.csv")

eigenvals(PCA.quimica)[1:2] #autovetor. Sempre tem que ser maior que o Broken-Stick
bstick(PCA.quimica)[1:2] #broken-stick
PCA.quimica$rotation[,1:2] #loadings. A import?nia de cada vari?vel para os eixos

PCAloadings2 <- data.frame(Variables = rownames(PCA.quimica$rotation), PCA.quimica$rotation)
SelLoad2<-PCAloadings2
SelLoad2[abs(SelLoad2[,3])>0.3,]
SelLoad2[,1][abs(SelLoad2[,2])>0.3]
SelLoad2[,3][abs(SelLoad2[,3])>0.3]
SelLoad2[,1][abs(SelLoad2[,3])>0.3]

unique(c(SelLoad2[,1][abs(SelLoad2[,2])>0.3], SelLoad2[,1][abs(SelLoad2[,3])>0.3]))->SelLoads2
SelLoad2[grepl(paste0(SelLoads2,collapse = "|"),SelLoad2$Variables),]->SelLoad2

cor.test(env$IIH_total,PCA.quimica$x[,1])
cor.test(env$IIH_total,PCA.quimica$x[,2])
cor(PCA.local$x[,1],PCA.quimica$x[,2])
summary(lm(PCA.local$x[,1]~PCA.quimica$x[,1]))

library(betareg)
data.frame(scale(quimicas),HII=env$IIH_total)->LMDATA
glm(HII~.,data=LMDATA,family="gaussian")->LM
step(LM)->LM
write.csv(LM$anova,"results/table s2.csv")
betareg(LM$formula,data=LMDATA)->LM
summary(LM)
write.csv(summary(LM)[[1]]$mean,"results/table s8.csv")

P1<-ggplot(quimicas, aes(x=Phos, y=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="glm",method.args = list(family = 'quasipoisson'),color="black")+
  theme_bw()+
  #scale_x_continuous(limits = c(min(VarNonCor$Nitg),max(VarNonCor$Nitg)+0.05))+
  labs(x="Total phosphorus (mg-P/L)", y="Habitat integrity index")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); P1

P2<-ggplot(quimicas, aes(x=Amm, y=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="glm",method.args = list(family = 'quasipoisson'),color="black")+
  theme_bw()+
  scale_x_continuous(expand = c(0,0.03))+
  labs(x="Ammonia (mg/L)", y="Habitat integrity index")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); P2

P3<-ggplot(quimicas, aes(x=Alum, y=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="glm",method.args = list(family = 'quasipoisson'),color="black")+
  theme_bw()+
  #scale_x_continuous(limits = c(min(VarNonCor$Nitg),max(VarNonCor$Nitg)+0.05))+
  labs(x="Aluminum (mg/L)", y="Habitat integrity index")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); P3

library(ggpubr)
tiff("raw_figures/Figure 3.tiff", units="in", width=7.5, height=5.5, res=600)
ggarrange(P1,P2,P3,nrow=2,ncol=2,labels = c("(a)","(b)","(c)"),label.x = 0.905, label.y = 0.98, common.legend = T,legend = "right",hjust = 0)
dev.off()

SelLoad2$Variables=c("BOD","Total\nPhosphorus","Nitrate","Ammonia","DO","Aluminum","Iron")

pca2<-ggplot(data.frame(PCA.quimica$x), aes(x = PC1, y = PC2)) +
  geom_hline(yintercept=0, linetype="dashed", alpha=0.5)+
  geom_vline(xintercept=0, linetype="dashed", alpha=0.5)+
  geom_point(aes(x = PC1, y = PC2, color=HII$IIH_total), size=3)+
  geom_segment(data = SelLoad2, aes(x = 0, y = 0, xend = (PC1*2.5),
                                    yend = (PC2*2.5)), arrow = arrow(length = unit(1/4, "picas")),color = "black", size=0.6,alpha=1) +
  geom_text_repel(data=SelLoad2,aes(x = (SelLoad2$PC1*2.6), y = (SelLoad2$PC2*2.5)),
                  label = SelLoad2$Variables, color="black", fontface="italic", size=3,alpha=1)+
  scale_x_continuous(name="PC1 (34.82%)")+
  scale_y_continuous(name="PC2 (20.95%)")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  theme(legend.text.align = 0,legend.title.align = 0,legend.key.size = unit(0.5,"cm"),legend.position ="right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(family = "sans",size=12, face="bold"), legend.text = element_text(family = "sans",size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); pca2

tiff("raw_figures/Figure 2.tiff", units="in", width=8, height=3, res=600)
ggpubr::ggarrange(pca1,pca2,nrow=1,ncol=2,labels = c("(a)","(b)"),label.x = 0.9, label.y = 0.98, common.legend = T,legend = "right",hjust = 0)
dev.off()

####### Três facetas da diversidade #####
library(vegan)
rowSums(decostand(comm[,1:ncol(comm)-1],"pa"))->rich
S <- apply(decostand(comm[,1:ncol(comm)-1],"pa")>0,1,sum)
diversity(comm[,1:ncol(comm)-1], index="shannon")/log(S)->even
library(adespatial)
beta.div(decostand(comm[,1:ncol(comm)-1],"pa"))->beta
beta$LCBD->dist
#rowSums(comm[,1:ncol(comm)-1])->Dens

cor(cbind(rich,even,dist))

#OS glm#
data.frame(rich=rich,scale(cbind(quimicas,HII=HII$IIH_total)))->LMData1
data.frame(rich=rich,(cbind(quimicas,HII=HII$IIH_total)))->NoPad1
glm(rich~.,family = "poisson",data=LMData1)->LM1
step(LM1,direction = "both")->LM1
write.csv(LM1$anova,"results/table s3.csv")
summary(LM1)
write.csv(summary(LM1)[[12]],"results/table s10.csv")
with(summary(LM1), 1 - deviance/null.deviance)
lm(LM1$formula,data=NoPad1)->LMNoPad1
summary(LMNoPad1)

data.frame(even=even,scale(cbind(quimicas,HII=HII$IIH_total)))->LMData2
data.frame(even=even,(cbind(quimicas,HII=HII$IIH_total)))->NoPad2
glm(even~.,family = "gaussian",data=LMData2)->LM2
step(LM2,direction = "both")->LM2
write.csv(LM2$anova,"results/table s4.csv")
betareg(LM2$formula,data=LMData2)->LM2_final
summary(LM2_final)
write.csv(summary(LM2_final)[[1]]$mean,"results/table s11.csv")
lm(LM2$formula,data=NoPad2)->LMNoPad2
summary(LMNoPad2)

data.frame(dist=dist,scale(cbind(quimicas,HII=HII$IIH_total)))->LMData3
glm(dist~.,family = "gaussian",data=LMData3)->LM3
step(LM3,direction = "both")->LM3
write.csv(LM3$anova,"results/table s5.csv")
betareg(LM3$formula,data=LMData3)->LM3_final
summary(LM3_final)
write.csv(summary(LM3_final)[[1]]$mean,"results/table s12.csv")

data.frame(Variables=c(rownames(summary(LM1)[[12]]),rownames(summary(LM2_final)[[1]]$mean),rownames(summary(LM3_final)[[1]]$mean)),rbind(summary(LM1)[[12]],summary(LM2_final)[[1]]$mean, summary(LM3_final)[[1]]$mean),index=c(rep("Richness",nrow(summary(LM1)[[12]])),rep("Evenness",nrow(summary(LM2_final)[[1]]$mean)),rep("Divergence",nrow(summary(LM3_final)[[1]]$mean))))->indexPlot
indexPlot[!indexPlot$Variables=="(Intercept)",]->indexPlot2

indexPlot2$lower_Estimate<-indexPlot2$Estimate-(indexPlot2[,3]*2)
indexPlot2$upper_Estimate<-indexPlot2$Estimate+(indexPlot2[,3]*2)

tiff("raw_figures/Figure 5.tiff", units="in", width=6.5, height=3, res=600)
ggplot(indexPlot2, aes(x = Variables, y = Estimate)) +
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.75)+
  geom_pointrange(aes(ymin = lower_Estimate, ymax = upper_Estimate), linewidth=0.8,size=0.5)+
  facet_wrap(.~factor(index,levels = c("Richness","Evenness","Divergence")))+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(name="Effect size")+
  scale_x_discrete(name="Explanatory variables",labels=c("Aluminum","Ammonia","Habitat integrity","Nitrate","Dissolved oxygen","Total phosphorus"))+
  theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))
dev.off()

Nitg_plot<-ggplot(quimicas, aes(x=Nitg, y=rich, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="glm",method.args = list(family = 'poisson'),color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  scale_x_continuous(limits = c(min(VarNonCor$Nitg),max(VarNonCor$Nitg)+0.05))+
  labs(x="Nitrate (mg/L)", y="Taxonomic richness")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));Nitg_plot

HII_plot<-ggplot(HII, aes(x=IIH_total, y=rich, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="glm",method.args = list(family = 'poisson'),color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  labs(x="Habitat integrity index", y="Taxonomic richness")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));HII_plot

OD_plot<-ggplot(quimicas, aes(x=OD, y=even, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="glm",method.args = list(family = 'quasipoisson'),color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  labs(x="Dissolved oxygen (mg/L)", y="Taxonomic evenness")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));OD_plot

tiff("raw_figures/Figure 3.tiff", units="in", width=4.4*2, height=3.3*1.75, res=600)
ggpubr::ggarrange(Nitg_plot,HII_plot,OD_plot,nrow=2,ncol=2,labels = c("(a)","(b)","(c)"),label.x = 0.915, label.y = 0.98, common.legend = T,legend = "right",hjust = 0)
dev.off()

library(reshape2)
data.frame(index=rep(dist,5),melt(data.frame("Total phosphorus"=quimicas$Phos,"Nitrate"=quimicas$Nitg,"Ammonia"=quimicas$Amm,"Dissolved oxygen"=quimicas$OD,"Aluminum"=quimicas$Alum)))->dist_plot
gsub("Dissolved.oxygen","Dissolved oxygen",dist_plot$variable)->dist_plot$variable

dist_plots<-lapply(unique(dist_plot$variable),function(x){
  ggplot(dist_plot[dist_plot$variable==x,], aes(x=value, y=index, color=rep(HII$IIH_total,1))) + 
    geom_point(size=2.5)+
    geom_smooth(method="glm",method.args = list(family = 'quasipoisson'),color="black")+
    # facet_wrap(.~factor(variable,levels=c("Aluminum","Ammonia","Dissolved oxygen","Phosphorus","Nitrate")),scale="free_x")+
    scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
    theme_bw()+
    labs(x=paste0(x," (mg/L)"), y="Taxonomic divergence")+ #Mudar o t?tulo da legenda
    theme(strip.text = element_text(size=12,family="sans",face = "bold"),strip.background = element_rect(fill="white"),legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))
  
})

tiff("raw_figures/Figure 4.tiff", units="in", width=13, height=6, res=600)
ggpubr::ggarrange(dist_plots[[1]],dist_plots[[3]],dist_plots[[2]],dist_plots[[4]],dist_plots[[5]],nrow=2,ncol=3,labels = c("(a)","(b)","(c)","(d)","(e)"),label.x = 0.915, label.y = 0.98, common.legend = T,legend = "right",hjust = 0)
dev.off()
