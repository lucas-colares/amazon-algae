# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############# Statistical analysis for Palheta et al. (2023) ##############
# # # # # # # # # # By: Lucas Colares & Leandra Palheta # # # # # # # # # #
# # # # # # # # # # # # # # # # 20-08-2023 # # # # # # # ## # # # # # # # #

read.csv("datasets/comm.csv",h=T,sep=";",row.names = 1)->comm
read.csv("datasets/env.csv",h=T,sep=";",row.names = 1)->env
read.csv("datasets/HII.csv",h=T,sep=";",row.names = 1)->HII

###### Environment ----
#### Linear models for environmental variables
quimicas=env[,1:9]
fisicas=HII[,1:(ncol(HII)-2)]
landscape=env[,11:15]
# Chemical block
write.csv(cor(scale(data.frame(quimicas,landscape))),"results/correlation matrix - chemical variables.csv")
VarNonCor<-data.frame(quimicas[,c("BOD","Phos","Nitg","Amm","Alum","Iron")],landscape[,c("floresta","solo_exposto","urbano")])
LMData1<-data.frame(scale(VarNonCor),IIH_total=HII$IIH_total)
LM1=lm(IIH_total~.,data=LMData1)
LM1=step(LM1,direction = "both")
summary(LM1)
write.csv(cheLM$anova,"results/stepwise matrix - chemical variables.csv")
write.csv(summary(cheLM)[[4]],"results/coefs matrix - chemical variables.csv")

prcomp(fisicas)->PCA.local
summary(PCA.local)
eigenvals(PCA.local)[1:2] #autovetor. Sempre tem que se rmaior que o Broken-Stick.
bstick(PCA.local)[1:2] #broken-stick
PCA.local$rotation[,1:2] #loadings. A import?nia de cada vari?vel para os eixos.

PCAloadings1 <- data.frame(Variables = rownames(PCA.local$rotation), PCA.local$rotation)
SelLoad1<-PCAloadings1
SelLoad1[,2][abs(SelLoad1[,2])>0.3]
SelLoad1[,1][abs(SelLoad1[,2])>0.3]
SelLoad1[,3][abs(SelLoad1[,3])>0.3]
SelLoad1[,1][abs(SelLoad1[,3])>0.3]

SelLoad1$Variables=c("LandInt","WidtInt","ComplRip","VegInt","RetInt","SedInt","BnkInt","BnkUnd","BotInt","FlwComp","AquVegInt","DetInt")

pca1<-ggplot(data.frame(PCA.local$x), aes(x = PC1, y = PC2)) +
  geom_hline(yintercept=0, linetype="dashed", alpha=0.8)+
  geom_vline(xintercept=0, linetype="dashed", alpha=0.8)+
  geom_point(aes(x = PC1, y = PC2, color=HII$IIH_total), size=3)+
  geom_segment(data = SelLoad1, aes(x = 0, y = 0, xend = (PC1*2),
                                    yend = (PC2*2)), arrow = arrow(length = unit(1/4, "picas")),color = "grey55", size=0.8,alpha=0.8) +
  geom_text_repel(data=SelLoad1,aes(x = (SelLoad1$PC1*2), y = (SelLoad1$PC2*2)),
                  label = SelLoad1$Variables, color="black", fontface="italic", size=3,alpha=1)+
  scale_x_continuous(name="PC1 (62.66%)")+
  scale_y_continuous(name="PC2 (21.09%)")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  theme(legend.text.align = 0,legend.title.align = 0,legend.key.size = unit(0.5,"cm"),legend.position ="right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(family = "sans",size=12, face="bold"), legend.text = element_text(family = "sans",size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); pca1

Urban_env_plot<-ggplot(VarNonCor, aes(x=urbano, y=HII$IIH_total, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="lm",color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  labs(x="Urban cover (%)", y="Habitat integrity index")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));Urban_env_plot

Nitg_env_plot<-ggplot(VarNonCor, aes(x=Nitg, y=HII$IIH_total, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="lm",color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  scale_x_continuous(limits = c(min(VarNonCor$Nitg),max(VarNonCor$Nitg)+0.05))+
  labs(x="Nitrate (mg/L)", y="Habitat integrity index")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));Nitg_env_plot

tiff("raw_figures/Figure 2_copy.tiff", units="in", width=4.5*3, height=3.3, res=300)
ggpubr::ggarrange(pca1,Nitg_env_plot,Urban_env_plot,nrow=1,ncol=3,labels = c("(a)","(b)","(c)"),label.x = 0.92, label.y = 0.985, align = "v",common.legend = T,legend = "right",hjust = 0)
dev.off()

####### Riqueza e densidade #####
rowSums(decostand(comm[,1:ncol(comm)-1],"pa"))->S
rowSums(comm[,1:ncol(comm)-1])->Dens

#OS glm#
glm(S~.,family = "poisson",data=LMData1)->LM2
step(LM2,direction = "both")->LM2
summary(LM2)

Nitg_plot<-ggplot(VarNonCor, aes(x=Nitg, y=S, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="lm",color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  scale_x_continuous(limits = c(min(VarNonCor$Nitg),max(VarNonCor$Nitg)+0.05))+
  labs(x="Nitrate (mg/L)", y="Number of species")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));Nitg_plot

Urban_plot<-ggplot(VarNonCor, aes(x=urbano, y=S, color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="lm",color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity\ngradient")+
  theme_bw()+
  labs(x="Urban cover (%)", y="Number of species")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));Urban_plot

lm(log10(Dens)~.,data=LMData1)->LM3
step(LM3,direction = "both")->LM3
summary(LM3)

P_plot<-ggplot(VarNonCor, aes(x=Phos, y=log10(Dens), color=HII$IIH_total)) + 
  geom_point(size=2.5)+
  geom_smooth(method="lm",color="black")+
  scale_color_gradient(low = "red",high = "green",name="Integrity gradient")+
  theme_bw()+
  scale_y_continuous(breaks=c(4,5,6),labels = c(10^4,10^5,10^6))+
  labs(x="Phosphorus (mg/L)", y="Density (ind/cm2)")+ #Mudar o t?tulo da legenda
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"));P_plot

tiff("raw_figures/Figure 3.tiff", units="in", width=4.5*3, height=3.3, res=900)
ggpubr::ggarrange(Nitg_plot,Urban_plot,P_plot,nrow=1,ncol=3,labels = c("(a)","(b)","(c)"),label.x = 0.92, label.y = 0.985, common.legend = T,legend = "right",hjust = 0)
dev.off()