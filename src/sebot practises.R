library(phytools);library (ape);library(nlme);library(diversitree)
library(geiger); library(OUwie);library(tidyverse)
# data exploration
# from curso sebot pcm S1 Contrastes independientes #####
read.csv("data/species.csv")%>%
  select(Taxon, PERoil, seedmass, ratio, ecology)%>%
  group_by(Taxon, ecology)%>%
  summarise(oilcontent = mean(PERoil), 
            seedmass = mean(seedmass),
            ratio = mean(ratio))%>%
  mutate(Taxon= gsub(" ","_", Taxon))%>%
  column_to_rownames(var="Taxon") -> oil

## create scatterplot
plot(oilcontent~seedmass,data=oil,xlab="Seed mass (mg)",ylab="Oil content (%)",
     pch=21,bg="gray",cex=1.2,log="xy",las=1,cex.axis=0.7,cex.lab=0.9,bty="n")

plot(ratio~seedmass,data=oil,xlab="Seed mass (mg)",ylab="Oil content (%)",
     pch=21,bg="gray",cex=1.2,log="xy",las=1,cex.axis=0.7,cex.lab=0.9,bty="n")

#modelo regresion estandard
fit.ols<-lm(log(oilcontent)~log(seedmass),data=oil)
fit.ols
summary(fit.ols) # contenido en aceites depende de la masa de la semilla

fit.ols2<-lm(ratio~log(seedmass),data=oil)
fit.ols2
summary(fit.ols2) # ratio insaturados saturados no depende de la masa de la semilla

fit.ols3<-lm(ratio~log(oilcontent),data=oil)
fit.ols3
summary(fit.ols3) # ratio de UFa/SFA muy dependiente del oil content

## create scatterplot
plot(oilcontent~seedmass,data=oil,xlab="Seed mass (mg)",ylab="Oil content (%)",
     pch=21,bg="gray",cex=1.2,log="xy",las=1,cex.axis=0.7,cex.lab=0.9,bty="n")
## add the line of best fit from lm
lines(oil$seedmass,exp(predict(fit.ols)),lwd=2,col="darkgray")

plot(ratio~seedmass,data=oil,xlab="Seed mass (mg)",ylab="Ratip UFA/SFA",
     pch=21,bg="gray",cex=1.2,log="xy",las=1,cex.axis=0.7,cex.lab=0.9,bty="n")
## add the line of best fit from lm
lines(oil$seedmass,exp(predict(fit.ols2)),lwd=2,col="darkgray")

# arbol filogenético 
tree<-read.tree("results/tree.tree")
## plot phylogeny of species
plotTree(tree,ftype="i",fsize=0.7,lwd=1) # 
## add node labels to the plotted tree
nodelabels(bg="white",cex=0.5,frame="circle") # 

## check to see if data and phylogeny match using
## geiger::name.check
chk<-name.check(tree,oil)
summary(chk)


## pull our home range and body mass as
## numeric vectors
oilcontent <-setNames(oil[,"oilcontent"],rownames(oil))
ratio <-setNames(oil[,"ratio"],rownames(oil))
seedmass<-setNames(oil[,"seedmass"],rownames(oil))

## compute PICs for home range and body size
pic.oilcontent<-pic(log(oilcontent),tree)
pic.ratio<-pic(ratio,tree)
pic.seedmass<-pic(log(seedmass),tree)

head(pic.oilcontent,n=20)
head(pic.seedmass,n=20)

## fit linear model to PICs without intercept
fit.pic<-lm(pic.oilcontent~pic.seedmass+0)
fit.pic
summary(fit.pic) #con contrastes independientes contenido de aceites no depende de la masa de la semilla

fit.pic2<-lm(pic.ratio~pic.seedmass+0)
fit.pic2
summary(fit.pic2) # ratios sigue sin depender de la masa de la semilla

fit.pic3<-lm(pic.ratio~pic.oilcontent+0)
fit.pic3
summary(fit.pic3) # ratios depende claramente del contenido total de aceites

## set margins
par(mar=c(5.1,5.1,1.1,1.1))
## graph scatterplot of contrasts
plot(pic.oilcontent~pic.seedmass,xlab="PICs for log(seed mass)",ylab="PICs for log(oil content)",pch=21,bg="gray",cex=1.2,las=1,cex.axis=0.7,cex.lab=0.9,bty="n")
## add gridlines to the plot
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
## reset graphing limits of the plot to the## x/y range of our PICs
clip(min(pic.seedmass),max(pic.seedmass),min(pic.oilcontent),max(pic.oilcontent))
## graph our fitted line
abline(fit.pic,lwd=2,col="darkgray")

## set margins
par(mar=c(5.1,5.1,1.1,1.1))
## graph scatterplot of contrasts
plot(pic.ratio~pic.oilcontent,xlab="PICs for log(oil content)",ylab="PICs for ratio",pch=21,bg="gray",cex=1.2,las=1,cex.axis=0.7,cex.lab=0.9,bty="n")
## add gridlines to the plot
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
## reset graphing limits of the plot to the## x/y range of our PICs
clip(min(pic.oilcontent),max(pic.oilcontent),min(pic.ratio),max(pic.ratio))
## graph our fitted line
abline(fit.pic3,lwd=2,col="darkgray")

# from curso sebot pcm S2 GLS #####

spp<-rownames(oil)
corBM<-corBrownian(phy=tree,form=~spp)
corBM

pgls.oil<-gls(log(oilcontent)~log(seedmass),data=oil,correlation=corBM)
summary(pgls.oil) # sontenido de aceites no depende de la masa de la semilla (igual que con los PIC)
pgls.ratio<-gls(log(ratio)~log(seedmass),data=oil,correlation=corBM)
summary(pgls.ratio) 

coef(fit.pic)
coef(pgls.oil)
abs(coef(fit.pic)[1] -coef(pgls.oil)[2])

corLambda <-corPagel(value=1,phy=tree,form=~spp) # la función calcula el "mejor" valor de lamba, el valor que maximiza la probabilidad de que nuestros datos a aparezcan
corLambda # lambda = 1, significa que la variación observada proviene de la filogenia (sigue modelo browniano)
#  Concluimos que existe una correlación evolutiva entre los dos rasgos.

pgls.Lambda<-gls(log(oilcontent)~log(seedmass), data=oil,correlation=corLambda)
summary(pgls.Lambda)
# ahora el contenido de aceite depende de la massa de la semilla (marginalmente significativo)
pgls.Lambda2<-gls(ratio~log(seedmass), data=oil,correlation=corLambda)
summary(pgls.Lambda2) # sigue sin ser significativo 
# probamos UFA/SFA ~ oil content
pgls.Lambda3<-gls(ratio~log(oilcontent), data=oil,correlation=corLambda)
summary(pgls.Lambda3) # el ratio depende del contenido de aceites

#ANOVA FILOGENÉTICO ##
oil.ancova<-gls(log(oilcontent)~log(seedmass)+ecology,data=oil,correlation=corBM)
anova(oil.ancova) # ecology marginalmente significativo

ratio.ancova<-gls(ratio~log(oilcontent)+ecology,data=oil,correlation=corBM)
anova(ratio.ancova) # ecology NO significativo, oil content si

# from curso sebot pcm S3 modelo browniano y UO ####
# Necesitamos ajustar dos parámetros. Por una lado el parámetro de varianza instantanea del modelo browniano, 
# que hemos llamada sigma2. Y en segundo lugar necesitamos también estimar el valor de nuestro fenotipo 
# bajo estudio en la raíz de la filogenia.

par(mfrow=c(1,1))

# ajuste entre los nombres de la filogenia y set de datos
name.check(tree, oil)

#We create a vector for the function fitContinuous
oilcontent<-oil[,"oilcontent"]
oilcontent

ratio <- oil[, "ratio"]
ratio

seedmass <- oil[, "seedmass"]
seedmass
#We asign names to phenotypes in the vector
names(oilcontent)<-rownames(oil)
head(oilcontent)

names(ratio)<-rownames(oil)
head(ratio)

names(seedmass)<-rownames(oil)
head(seedmass)

## histogram  on original scale
hist(oilcontent,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(oilcontent),max(oilcontent),length.out=12))
mtext("(a)",adj=0,line=1)
mtext("rate",side=1,line=4,cex=0.9)

hist(ratio,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(ratio),max(ratio),length.out=12))
mtext("(a)",adj=0,line=1)
mtext("rate",side=1,line=4,cex=0.9)

hist(seedmass,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(seedmass),max(seedmass),length.out=12))
mtext("(a)",adj=0,line=1)
mtext("rate",side=1,line=4,cex=0.9)

## histogram of mutation accumulation rates on log scale
log_oilcontent<-log(oilcontent)
hist(log_oilcontent,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(log_oilcontent),max(log_oilcontent),length.out=12))
mtext("(b)",adj=0,line=1)
mtext("ln(rate)",side=1,line=4,cex=0.9)

log_seedmass<-log(seedmass)
hist(log_seedmass,main="",las=2,xlab="",cex.axis=0.7,cex.lab=0.9,breaks=seq(min(log_seedmass),max(log_seedmass),length.out=12))
mtext("(b)",adj=0,line=1)
mtext("ln(rate)",side=1,line=4,cex=0.9)

## fit Brownian motion model using fitContinuous
fitBM_oil<-fitContinuous(tree,log_oilcontent)
fitBM_oil
# sigma = 0.011159 (por cada unidad de tiempo contenido de aciete cambia con una varianza de 0.649..); x0 (z0) = 1.888710
fitBM_ratio<-fitContinuous(tree,ratio)
fitBM_ratio
# sigma = 0.060163; x0 (z0) = 6.777468
fitBM_seedmass<-fitContinuous(tree,log_seedmass)
fitBM_seedmass
# sigma = 0.0205; x0 (z0) = 3.134499

## compute phylogenetic signal, lambda, 
phylosig(tree,log_oilcontent,method="lambda") #lambda = 0.506108 logL -35.7
phylosig(tree,ratio,method="lambda") # lambda = 0.795678 logL=-70.48
phylosig(tree,log_seedmass,method="lambda") # lamba = 0.988725 logL = -53.2149

## test for significant phylogenetic signal, lambda,
## in each of our traits
lambda_oil<-phylosig(tree,log_oilcontent,method="lambda",test=TRUE)
lambda_oil # lambda significantly different to 0 
lambda_ratio<-phylosig(tree, ratio,method="lambda",test=TRUE)
lambda_ratio # lambda significantly different to 0 
lambda_seedmass<-phylosig(tree,log_seedmass,method="lambda",test=TRUE)
lambda_seedmass # lambda significantly different to 0 

## our plot area into 1 column and two rows
par(mfrow=c(3,1),mar=c(5.1,4.1,2.1,2.1),cex=0.8)
## plot the likelihood surfaces of lambda for each of our
## two traits
plot(lambda_oil,las=1,cex.axis=0.9,bty="n",xlim=c(0,1.1))
mtext("(a) oil content",line=1,adj=0)
plot(lambda_ratio,las=1,cex.axis=0.9,bty="n",xlim=c(0,1.1))
mtext("(b) ratio UFA/SFA",line=1,adj=0)
plot(lambda_seedmass,las=1,cex.axis=0.9,bty="n",xlim=c(0,1.1))
mtext("(c) seed mass",line=1,adj=0)

#MODELO ORNSTEIN-UHLENBECK (OU)
#extensión e del movimiento browniano, con un parámetro adicional de selección (?? o alpha) que describe la tendencia a 
# evolucionar hacia un valor óptimo (??). Debido a que este modelo implica un cambio evolutivo hacia un valor particular, 
# se interpreta con mayor frecuencia como un modelo de evolución adaptativa en el que ?? (alpha) corresponde a la fuerza 
# de la selección natural y ?? a la posición del óptimo. 

## fit OU model 
fitOU_oil<-fitContinuous(tree,log_oilcontent,model="OU")
fitOU_oil
# alpha = 0.0245, sigma = 0.027795, x0 = 1.88454
fitOU_ratio<-fitContinuous(tree,ratio,model="OU")
fitOU_ratio
# alpha = 0.01305, sigma = 0.120529, x0 = 6.872
fitOU_seedmass<-fitContinuous(tree,log_seedmass,model="OU")
fitOU_seedmass
# alpha = 0.001447 sigma = 0.022766, x0 = 3.1326

# Comparación modelos browniano vs OU
## accumulate AIC scores from our three models into a vector
#oil content (OU better model) HAY un óptimo al que evolucionar
aic_oil<-setNames(c(AIC(fitBM_oil),AIC(fitOU_oil)),c("BM","OU"))
aic_oil # menos aic mejor
aic.w(aic_oil) #weighted = probabilidad del modelo
# ratio UFA/SFA (OU better model) HAY un óptimo al que evolucionar
aic_ratio<-setNames(c(AIC(fitBM_ratio),AIC(fitOU_ratio)),c("BM","OU"))
aic_ratio # menos aic mejor
aic.w(aic_ratio)#weighted = probabilidad del modelo
# seed mass (Browniano better model) No hay un óptimo al que evolucionar
aic_seedmass<-setNames(c(AIC(fitBM_seedmass),AIC(fitOU_seedmass)),c("BM","OU"))
aic_seedmass # menos aic mejor
aic.w(aic_seedmass) #weighted = probabilidad del modelo

# from curso sebot pcm s4 modelo browniano y UO con multiple tasas  ecology inserted in map! ####
## run phylogenetic PCA and print the results
# only numeric variables for PCA
par(mfrow=c(1,1)) 
dev.off()
oil%>%
  select(oilcontent, seedmass, ratio)->oil.pca
phylo.pca<-phyl.pca(tree,oil.pca)
plot(phylo.pca)
print(phylo.pca)

## create our OUwie data frame
ouwie.data<-data.frame(Genus_species=rownames(scores(phylo.pca)),Reg=oil[rownames(scores(phylo.pca)),],X=as.numeric(scores(phylo.pca)[,2]))
head(ouwie.data,n=10)

mtree<-make.simmap(tree,ecology,model="ARD")
cols<-setNames(rainbow(n=2),levels(species.ecology[,1]))
plot(mtree,cols,lwd=2,ftype="i",fsize=0.8,ylim=c(-4,82),outline=T)
add.simmap.legend(colors=cols,prompt=FALSE,x=0,y=-2,vertical=FALSE,fsize=0.9)

## fit standard, one-rate Brownian model specified by model statements
fitBM<-OUwie(mtree,ouwie.data,model="BM1",simmap.tree=T)
fitBM

# fit multi-rate Brownian model specified by model statements
fitBMS<-OUwie(mtree,ouwie.data,model="BMS",simmap.tree=T,root.station=FALSE)
fitBMS

# fit multi-regime OU model
fitOUM<-OUwie(mtree,ouwie.data,model="OUM",simmap.tree=T,root.station=FALSE)
fitOUM

## extracting AIC scores
aic<-setNames(c(fitBM$AIC,fitBMS$AIC,fitOUM$AIC),c("BM1","BMS","OUM"))
aic
## compute Akaike weights
aic.w(aic)
# El modelo OU con múltiples optimos o Multi rate brownian model son mejores
ecology<-setNames(oil.data[,1],rownames(oil.data))
oil%>%
  select(ecology)%>%
  mutate(ecology = as.factor(ecology)) -> species.ecology
## get the tip state for  for each species
tips<-getStates(mtree,"tips")
## set these tip states to have the colors using the
## color scheme 
tip.cols<-cols[tips]
## plot tree with adjacent barplot using
## phytools::plotTree.barplot using PCA2 
plotTree.barplot(mtree,scores(phylo.pca)[,2],args.plotTree=list(fsize=0.85),args.barplot=list(col=tip.cols,xlab=expression(paste("PC2 (",""%up%"","oil content, ",""%down%"","lamellae)",sep="")),cex.lab=0.8))
## add an informative legend
legend("topright",levels(species.ecology[,1]),pch=22,pt.bg=cols,pt.cex=1.5,cex=0.9)

# Use trait values OIL CONTENT
plotTree.barplot(mtree,oilcontent,args.plotTree=list(fsize=0.85),args.barplot=list(col=tip.cols,xlab="Oil content (%)",cex.lab=0.8))
## add an informative legend
legend("topright",levels(species.ecology[,1]),pch=22,pt.bg=cols,pt.cex=1.5,cex=0.9)

# Use trait values ratio
plotTree.barplot(mtree,ratio,args.plotTree=list(fsize=0.85),args.barplot=list(col=tip.cols,xlab="Ratio UFA/SFA",cex.lab=0.8))
## add an informative legend
legend("right",levels(species.ecology[,1]),pch=22,pt.bg=cols,pt.cex=1.5,cex=0.9)

# Use trait values seed mass
plotTree.barplot(mtree,seedmass,args.plotTree=list(fsize=0.85),args.barplot=list(col=tip.cols,xlab="Seed mass (mg)",cex.lab=0.8))
## add an informative legend
legend("topright",levels(species.ecology[,1]),pch=22,pt.bg=cols,pt.cex=1.5,cex=0.9)

# modelos brownianos ocn multiples variables y correlaciones evolutivas
# H: ser especialista o generalista determina el contenido en aceites, ratio o seed mass
## convert numerical trait data to a matrix
oil.num<-as.matrix(oil.pca[, c("oilcontent", "ratio")])
## fit multi-regime multivariate Brownian model
fitMV<-evol.vcv(mtree,oil.num)
fitMV
# modelo más simple nos da 6 valores sigma (todas las combinaciones)
#log (L) = ajuste del modelo! más grande mejor! El modelo con dos matrics separadas para especialista y generalista mejor
#Interpretacion mediante correlaciones

cov2cor(fitMV$R.single)
cov2cor(fitMV$R.multiple[,,"Specialist"]) 
cov2cor(fitMV$R.multiple[,,"Generalist"]) # en generalistas mas correlacion (<0.9)entrer oil content y ratio

## modify the margins of the plot area and adjust axis
## and label font sizes
par(cex.lab=0.7,cex.axis=0.6,mar=c(5.1,4.1,1.1,2.1))
## plot the phylomorphospace without the mapped
## regimes
phylomorphospace(as.phylo(mtree),oil.num,ftype="off",lwd=4,xlab="oil content",ylab="ratio",node.size=c(0,0),bty="n")
## add the phylomorphospace projection with the
## mapped regimes
phylomorphospace(mtree,oil.num,ftype="off",colors=cols,lwd=2,node.by.map=TRUE,xlab="oil content",ylab="ratio",node.size=c(0,1.3),add=TRUE)
## add a legend
legend("topleft",c("Generalists","Specialists"),pch=22,pt.bg=cols,pt.cex=1.5,cex=0.7)

# HETEROGENEIDAD EVOLUTIVA DEL MODELO BROWNIANO 
## fit single-rate model (no rate shift) modelo browniano 1 sigma
fit1<-rateshift(tree,log_oilcontent) #log_oilcontent   ratio log_seedmass
## fit two-rate model (one rate shift) modelo browniano 2 sigma
fit2<-rateshift(tree,log_oilcontent,nrates=2)
## fit three-rate model (two rate shifts) modelo browniano 3 sigma
fit3<-rateshift(tree,log_oilcontent,nrates=3)
fit3
## fit EB model using geiger::fitContinuous tasa de sigma desciende con el tiempo
fitEB<-fitContinuous(tree,log_oilcontent,model="EB")
## compile our results into a list, sorted by
## the number of parameters estimated
fits<-list(fit1,fitEB,fit2,fit3) #
## create a table summarizing model fits
data.frame(model=c("BM","EB","two-rate","three-rate"),logL=sapply(fits,logLik),k=sapply(fits,function(x) attr(logLik(x),"df")),AIC=sapply(fits,AIC),weight=unclass(aic.w(sapply(fits,AIC))))

# for log_oil content 2 rates model is better (3 rates model gives an error)
# for ratio 3 rates model is better
# for seed mass BM is the better model (corcordingly bc it has a lambda of 0.98)

### visualization models!
## compute the total height of our  tree
h<-max(nodeHeights(tree))
## split our plot window into eight panels
par(mfrow=c(4,2))
## panel a) single-rate model graphed on the tree
plot(fit1,mar=c(1.1,4.1,2.1,0.1),ftype="i",fsize=0.5,col="gray")
mtext("(a)",adj=0,line=0)
## panel b) line graph of single-rate model
par(mar=c(4.1,4.1,2.1,1.1))
plot(NA,xlim=c(0,h),ylim=c(0,12),xlab="time",ylab=expression(sigma^2),bty="n")
lines(c(0,h),rep(fit1$sig2,2),lwd=3,col="gray")
mtext("(b)",adj=0,line=0)
## panel c) compute EB model and graph it on the tree
## calculate sigma^2 through time under fitted model
s2<-fitEB$opt$sigsq*exp(fitEB$opt$a*seq(h/200,h-h/200,length.out=100))
s2.index<-round((s2-min(s2))/diff(range(s2))*100)+1
## use make.era.map to paint fitted EB model onto tree
tmp<-make.era.map(tree,setNames(seq(0,h,length.out=101),s2.index))
## set colors for graphing
cols<-setNames(gray.colors(101,0.9,0),1:101)
## plot tree
plot(tmp,cols,mar=c(1.1,4.1,2.1,0.1),ftype="i",ylim=c(-0.1*Ntip(tree),Ntip(tree)),fsize=0.5)
## add color bar legend
add.color.bar(leg=0.5*h,cols=cols,prompt=FALSE,x=0,y=-0.05*Ntip(tree),lims=round(range(s2),3),title=expression(sigma^2))
mtext("(c)",adj=0,line=0)
## panel d) line graph of EB model
par(mar=c(4.1,4.1,2.1,1.1))
plot(NA,xlim=c(0,h),ylim=c(0,12),xlab="time",ylab=expression(sigma^2),bty="n")
lines(seq(0,h,length.out=100),s2,lwd=3,col="gray")
mtext("(d)",adj=0,line=0)
## panel e) two-rate model projected on the tree
plot(fit2,mar=c(1.1,4.1,2.1,0.1),ftype="i",fsize=0.5,col=cols)
mtext("(e)",adj=0,line=0)
## panel f) line graph of two-rate model
par(mar=c(4.1,4.1,2.1,1.1))
plot(NA,xlim=c(0,h),ylim=c(0,12),xlab="time",ylab=expression(sigma^2),bty="n")
lines(c(0,fit2$shift,h),c(fit2$sig2,fit2$sig2[2]),type="s",lwd=3,col="gray")
mtext("(f)",adj=0,line=0)
## panel g) three-rate model projected on the tree
plot(fit3,mar=c(1.1,4.1,2.1,0.1),ftype="i",fsize=0.5,col=cols)
mtext("(g)",adj=0,line=0)
## panel h) line graph of three-rate model
par(mar=c(4.1,4.1,2.1,1.1))
plot(NA,xlim=c(0,h),ylim=c(0,12),xlab="time",ylab=expression(sigma^2),bty="n")
lines(c(0,fit3$shift,h),c(fit3$sig2,fit3$sig2[3]),type="s",lwd=3,col="gray")
mtext("(h)",adj=0,line=0)

####ESTADÍSTICA BAYESIANA PARA LA COMPUTACIÓN DE MODELOS COMPLEJOS (paquetes con problema de instalación)
## load devtools
require(devtools)
## install bayou from GitHub
install_github("uyedaj/bayou")

library(bayou)
# from curso sebot pcm s5 Modelo markov , traits discretos ####
# por el momento no aplicable a nuestro conjunto de datos!
# from curso sebot pcm S6 modelo Pagel con caracteres discretos ####
par(mfrow=c(1,1))
# caracteres discretos binarios
oil%>%
  select(ecology)%>%
  mutate(ecology = as.factor(ecology)) -> species.ecology

## plot the tree with adjacent data matrix
object<-plotTree.datamatrix(tree,species.ecology,fsize=0.5,yexp=1,header=FALSE,xexp=1.45,palettes=c("YlOrRd","PuBuGn"))
## add a legend for trait 1
leg<-legend(x="topright",names(object$colors$ecology),cex=0.7,pch=22,pt.bg=object$colors$ecology,pt.cex=1.5,bty="n",title="Ecology")

# from curso sebot pcm S7  reconstruccion estados ancestrales + mapeo estocastico ####
# caracteres continuos ##
## read tree from file
print(tree,printlen=2)
oil%>%
  mutate(ecology = as.factor(ecology)) -> oil.data
head(oil.data)

head(log_oilcontent)
head(ratio)
head (log_seedmass)

## estimate ancestral states using fastAnc
fit.log_oilcontent<-fastAnc(tree,log_oilcontent,vars=TRUE,CI=TRUE)
print(fit.log_oilcontent,printlen=10)

fit.ratio<-fastAnc(tree,ratio,vars=TRUE,CI=TRUE)
print(fit.ratio,printlen=10)

fit.log_seedmass<-fastAnc(tree,log_seedmass,vars=TRUE,CI=TRUE)
print(fit.log_seedmass,printlen=10)

## plot eel phylogeny using plotTree
plotTree(tree,ftype="i",fsize=0.5,lwd=1)
## add node labels for reference
labelnodes(1:tree$Nnode+Ntip(tree),1:tree$Nnode+Ntip(tree),interactive=FALSE,cex=0.5)

#oil content (log transformed)
## compute "contMap" object
oil.contMap<-contMap(tree,log_oilcontent,plot=T,lims=c(0.1,3))
## change the color gradient to a custom gradient
oil.contMap<-setMap(oil.contMap,c("white","orange","black"))
## plot "contMap" object
plot(oil.contMap,sig=2,fsize=c(0.8,0.9),lwd=c(4,5),leg.txt="log(oil content %)")
## add error bars
errorbar.contMap(oil.contMap,lwd=8)

#oil content 
oilcontent<-oil[,"oilcontent"]
oilcontent
names(oilcontent)<-rownames(oil)
head(oilcontent)

## compute "contMap" object
oil.contMap<-contMap(tree,oilcontent,plot=T,lims=c(0.1,35))
## change the color gradient to a custom gradient
oil.contMap<-setMap(oil.contMap,c("white","orange","black"))
## plot "contMap" object
plot(oil.contMap,sig=2,fsize=c(0.8,0.9),lwd=c(4,5),leg.txt="oil content %")
## add error bars
errorbar.contMap(oil.contMap,lwd=8)

# ratio 
## compute "contMap" object
ratio.contMap<-contMap(tree,ratio,plot=T,lims=c(2,12))
## change the color gradient to a custom gradient
#ratio.contMap<-setMap(ratio.contMap,c("turquoise","green","forestgreen"))
## plot "contMap" object
plot(ratio.contMap,sig=2,fsize=c(0.8,0.9),lwd=c(4,5),leg.txt="Ratio UFA/SFA")
## add error bars
errorbar.contMap(ratio.contMap,lwd=8)

# seed mass
## compute "contMap" object
log_seedmass.contMap<-contMap(tree,log_seedmass,plot=T,lims=c(0.1,6))
## change the color gradient to a custom gradient
log_seedmass.contMap<-setMap(log_seedmass.contMap,c("white","orange","black"))
## plot "contMap" object
plot(log_seedmass.contMap,sig=2,fsize=c(0.8,0.9),lwd=c(4,5),leg.txt="log (seed mass(mg)")
## add error bars
errorbar.contMap(log_seedmass.contMap,lwd=8)

## caracteres discretos ###

## extract ecology as a vector
ecology<-setNames(oil.data[,1],rownames(oil.data))
## set colors for plotting
cols<-setNames(c("red","lightblue"),levels(ecology))
## plot the tree & data
plotTree.datamatrix(tree,as.data.frame(ecology),colors=list(cols),header=FALSE,fsize=0.85)
## add legend
legend("topright",legend=levels(ecology),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)

## fit ER model
fitER<-fitMk(tree,ecology,model="ER")
## fit ARD model
fitARD<-fitMk(tree,ecology,model="ARD")
## fit specialist -> generalist model
fit01<-fitMk(tree,ecology,model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
## fit generalist -> specialist model
fit10<-fitMk(tree,ecology,model=matrix(c(0,0,1,0),2,2,byrow=TRUE))
## extract AIC values for each model
aic<-c(AIC(fitER),AIC(fitARD),AIC(fit01),AIC(fit10))
## print summary table
data.frame(model=c("ER","ARD","bite->suction","suction->bite"),logL=c(logLik(fitER),logLik(fitARD),logLik(fit01),logLik(fit10)),AIC=aic,delta.AIC=aic-min(aic))

# el modelo con mejor ajuste es el ARD (all rates different)

library(corHMM)
## create new data frame for corHMM
oil.data2<-data.frame(Genus_sp=names(ecology),ecology=as.numeric(ecology)-1)
head(oil.data2,n=10)
## estimate marginal ancestral states under a ER model
##warning this make take a while
fit.marginal<-corHMM(tree,oil.data2,node.states="marginal",rate.cat=1,model ="ARD", rate.mat=matrix(c(NA,1,1,NA),2,2))
fit.marginal
head(fit.marginal$states)

## plot the tree & data
plotTree.datamatrix(tree,as.data.frame(ecology),colors=list(cols),header=FALSE,fsize=0.85)
## add legend
legend("topright",legend=levels(ecology),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)
## add node labels showing marginal ancestral states
nodelabels(pie=fit.marginal$states,piecol=cols,cex=0.5)

# mucha incertidumbre 
## Mapeo estocásticos ###

## generate one stochastic character history
mtree<-make.simmap(tree,ecology,model="ARD")
## plot single stochastic map
plot(mtree,cols,fsize=0.85,ftype="i",lwd=2,offset=0.4,ylim=c(-1,Ntip(tree)))
## add legend
legend("bottomleft",legend=levels(ecology),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)

## generate 1,000 stochastic character maps in which
## the transition rate is sampled from its posterior
## distribution
##warning, this may take a while!!!!
mtrees<-make.simmap(tree,ecology,model="ARD",nsim=1000,Q="mcmc",vQ=0.01,prior=list(use.empirical=TRUE),samplefreq=10)
mtrees

# set plot margins
par(mar=c(5.1,4.1,2.1,2.1))
## create a plot of the posterior density from stochastic
## mapping
plot(d<-density(sapply(mtrees,function(x) x$Q[1,2]),bw=0.005),bty="n",main="",xlab="q",xlim=c(0,0.5),ylab="Posterior density from MCMC",las=1,cex.axis=0.8)
polygon(d,col=make.transparent("blue",0.25))
## add line indicating ML solution for the same parameter
abline(v=fit.marginal$solution[1,2])
text(x=fit.marginal$solution[1,2],y=max(d$y),"MLE(q)",pos=4)

## create a 10 x 10 grid of plot cells
par(mfrow=c(10,10))
## graph 100 stochastic map trees, sampled evenly from
## our set of 1,000
null<-sapply(mtrees[seq(10,1000,by=10)],plot,colors=cols,lwd=1,ftype="off")

par(mfrow=c(1,1))
## compute posterior probabilities at nodes
pd<-summary(mtrees)
pd

## create a plot showing PP at all nodes of the tree
plot(pd,colors=cols,fsize=0.85,ftype="i",lwd=2,offset=0.4,ylim=c(-1,Ntip(tree)),cex=c(0.5,0.3))
## add a legend
legend("bottomleft",legend=levels(ecology),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)

# comparison marginal probabily with posterior probability
## set margins
par(mar=c(5.1,4.1,2.1,2.1))
## graph marginal ancestral states and posterior
## probabilities from stochastic mapping
plot(fit.marginal$states[,1],pd$ace[1:tree$Nnode],pch=21,cex=1.2,bg="grey",xlab="Marginal scaled likelihoods",ylab="Posterior probabilities",bty="n",las=1,cex.axis=0.8)
lines(c(0,1),c(0,1),col="blue",lwd=2)

## create a "densityMap" object
oil.densityMap<-densityMap(mtrees,states=levels(ecology)[2:1],plot=FALSE)

## update color gradient
## this make take a while!!!!
oil.densityMap<-setMap(oil.densityMap,cols[2:1])
## plot it, adjusting the plotting parameters
plot(oil.densityMap,fsize=c(0.8,0.9),lwd=c(5,4))

# from curso sebot pcm s8 LINAJES a través del tiempo ####
# resultados no fiables porque se necesita un elevado muestreo de species dentro del clado
# no es el caso de nuestros datos, de momento no aplicable.
# from curso sebot pcm 9 Tasas de diversificacion que varían a lo largo del tiempo ####
# mismo comentario que el anterior punto, no tenemos datos completos
# from curso sebot pcm 10 Diversificación dependiente del estado del caracter
# phylogenetic tree
tree
# trait data
oil
#ecology data # 1 = specialist, 2 = generalist
oil.data2%>%
  column_to_rownames(var="Genus_sp")->hab
hab<-setNames(hab[,1],rownames(hab))
## plot our tree
plotTree(tree,ftype="i",fsize=0.85,offset=0.5)
## add tip labels
tiplabels(pie=to.matrix(hab,0:1)[tree$tip.label,],piecol=c("white","black"),cex=0.4)
## create legend
legend("bottomleft",c("Generalist","Specialist"),pch=21,pt.cex=1.6,cex=0.8,bty="n",pt.bg=c("white","black"))

## make BiSSE likelihood function
bisse.model<-make.bisse(tree,hab)
## find reasonable parameter values for
## optimization
p<-starting.point.bisse(tree)
p
## optimize BiSSE model
bisse.mle<-find.mle(bisse.model,p)
bisse.mle

## create constrained null model
bissenull.model<-constrain(bisse.model,lambda1~lambda0,mu1~mu0)
## optimize null model
bissenull.mle<-find.mle(bissenull.model,p[c(-2,-4)])

coef(bissenull.mle)
logLik(bissenull.mle)
## run likelihood-ratio test
bisseAnova<-anova(bisse.mle,null=bissenull.mle)
bisseAnova

aicw(setNames(bisseAnova$AIC,rownames(bisseAnova)))
# No diferencias entre tasas de diversificación especialistas/generalistas

#Ahora vamos a ajustar nuestro modelo BiSSE con MCMC bayesiano.
prior<-make.prior.exponential(1/(2*0.4))
prior

## run Bayesian MCMC
bisse.mcmc<-mcmc(bisse.model,bisse.mle$par,nsteps=1000,prior=prior,w=0.1,print.every=100)
bisse.mcmc
##subdivide plot and set margins
par(mfrow =c(1,2),mar=c(5.1,4.1,3.1,2.1))
##set colors for plotting
col <-setNames(c("red","blue"), c("Generalist", "Specialist"))
## create graph of posterior sample for lamdda
profiles.plot(bisse.mcmc[,c("lambda0","lambda1")],col.line=col,las=1,bty="n",xlab=expression(lambda),cex.axis=0.7)
## add legend & panel label
legend("topright",names(col),pch=15,col=col,pt.cex=1.5,bty="n",cex=0.7)
mtext("a)",line=0.5,adj=0)
## create graph of posterior sample for mu
profiles.plot(bisse.mcmc[,c("mu0","mu1")],col.line=col,las=1,bty="n",xlab=expression(mu),cex.axis=0.7)
## add legend & panel label
legend("topright",names(col),pch=15,col=col,pt.cex=1.5,bty="n",cex=0.7)
mtext("b)",line=0.5,adj=0)

### lambda
sum(bisse.mcmc$lambda1>bisse.mcmc$lambda0)/length(bisse.mcmc$lambda1)
# No hay diferencias en tasas de especiacion entre especialistas y generalistas
## mu
sum(bisse.mcmc$mu1>bisse.mcmc$mu0)/length(bisse.mcmc$mu1)
# 
# MODELO HISSE (too complicated for my data)
# MODELO QUASSE tasa de especiación y/o extinción varian en función un caracter cuantitativo continuo (Fitzjohn 2010)