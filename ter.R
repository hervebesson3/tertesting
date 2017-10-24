#regression linéaire simple puis multiple
#install.packages("nortest")
#install.packages("lmtest")
#install.packages("HH")
#install.packages("FactorMineR")
library(HH)
library(xlsx)
library(lmtest)
library(nortest)
library(FactoMineR)
library(ggplot2)
library(lattice)
library(PerformanceAnalytics)
library(bios2mds)
library(vegan)
library(bio3d)

#1er dataframe
#datavariable<-read.table("Cohorte vascu Marie R.txt",sep = "\t",header = TRUE)
datavariable<-read.xlsx("try.xlsx",sheetIndex=1, header=TRUE, colClasses="character")
#gene <- read.table("Cohorte vascu Marie R.xls",sep="~",header = T)
rownames(datavariable)<-datavariable[,1]

datavariable<-datavariable[,-1]
datavariable[1:11,]
datavariables<-t(datavariable)
attach(datavariable)

#pairs pour all
datavariable[, 3] <- sapply(datavariable[, 3], as.numeric)
datavariable[, 1] <- sapply(datavariable[, 1], as.numeric)
datavariable[, 11] <- sapply(datavariable[, 11], as.numeric)
mydata <- datavariable[, c(1,2,3,4,5,6,7,8,9,10,11)]
chart.Correlation(mydata, histogram=TRUE, pch=19)

#2nd data frame

#read.table("Rj.txt", sep="\t", header=TRUE)

tableau<-read.table("Rj.txt", sep="\t", header=TRUE)
tab<-tableau[,c(1:20)]

##Go through each row and determine if a value is zero
row_sub = apply(tab, 1, function(row) all(row !=0 ))
##Subset as usual
datagenes<-tab[row_sub,]
datagene<-datagenes[1:200,c(2:15,17:20)]
datageness<-t(datagene)
sumGB<-apply(datagenes,1,sum)
sumGB1<-t(sumGB)

#age
fit<-lm(datageness~age, data=datavariable)
slope1<-coef(fit)
slope1
#sexe
fit<-lm(datageness~sexe, data=datavariable)
slope2<-coef(fit)
slope2
#reads
fit<-lm(datageness~reads, data=datavariable)
slope3<-coef(fit)
slope3
#concentration
fit<-lm(datageness~concentration, data=datavariable)
slope4<-coef(fit)
slope4
#RIN
fit<-lm(datageness~RIN, data=datavariable)
slope5<-coef(fit)
slope5


impact<-rbind(slope1,slope2,slope3,slope4,slope5)
impact<-as.data.frame(impact[c(2,4,6,8,10),])
cmon<-colSums(impact)
cmon

intercept<-rbind(slope1,slope2,slope3,slope4,slope5)
intercept<-as.data.frame(intercept[c(1,3,5,7,9),])
cmon1<-colSums(intercept)
cmon1

summary(fit)
attach(fit)
tygh<-fit$coefficients[,]
tygh
dummy.coef(fit)
summary(fit)[[4]]



#PCA
data<-read.table("try.xlsx", sep="\t", header=TRUE, row.names=1, check.names=F)
library(FactoMineR)
pca<-PCA(datageness)
summary(pca)
str(summary(pca))
library(ade4)
library(factoextra)
dudi.pca(datageness, center = TRUE,  scale = TRUE, 
         scannf = TRUE, nf = 20)
res.pca <- dudi.pca(datageness, scannf = FALSE, nf = 20)
summary(res.pca)
eig.val <- get_eigenvalue(res.pca)
head(eig.val)

#pc correlation :

#To examine variability of all numeric variables
sapply(datageness[],var)
range(sapply(datageness[],var))
# maybe this range of variability is big in this context.
#Thus, we will use the correlation matrix
#For this, we must standardize our variables with scale() function:
iris.stand <- as.data.frame(scale(datageness[,]))
sapply(iris.stand,sd) 

#If we use prcomp() function, we indicate 'scale=TRUE' to use correlation matrix
acp1 <- prcomp(iris.stand,scale=T)
#it is just the same that: prcomp(iris[,1:4],scale=T) and prcomp(iris.stand)
#similar with princomp(): princomp(iris.stand, cor=T)
acp1
summary(acp1)
#This gives us the standard deviation of each component, and the proportion of variance explained by each component.
#The stand




#plot of variance of each PCA.
#It will be useful to decide how many principal components should be retained.
screeplot(acp1, type="lines",col=3)

acp1$rotation

summary(acp1)
names(summary(acp1))



head(res.pca$co)

s.corcircle(res.pca$co)
fviz_pca_var(res.pca)
fviz_pca_contrib(res.pca, choice = "var", axes = 1)
fviz_pca_contrib(res.pca, choice = "var", axes = 2)
fviz_pca_contrib(res.pca, choice = "var", axes = 5)
head(res.pca$li)

prcomp(res.pca$co)

head(res.pca$co)
pc<-head(res.pca$co)[,1:5]
pc
prcomp(acp1$rotation)
holein1<-merge(pc, datavariable)

plot(acp1$rotation[,2], acp1$rotation[,3], pch = 19,  
     xlab="PC1 - 41.2%",ylab="PC2 - 18.4%")

abline(h=0, v=0, lty = 2)

text(acp1$rotation[,1], acp1$rotation[,2], labels=rownames(acp1$rotation),
     cex=0.7, pos = 3)

#plot of variance of each PCA.
#It will be useful to decide how many principal components should be retained.
screeplot(acp1, type="lines",col=3)

acp1$rotation

summary(acp1)
names(summary(acp1))


pc<-head(pca)
head(acp2)

prcomp(acp1$rotation)
holein1<-merge(acp2, datavariable, by=0, all = T)

help(pca)

plot(acp1$rotation[,2], acp1$rotation[,3], pch = 19,  
     xlab="PC1 - 41.2%",ylab="PC2 - 18.4%")

abline(h=0, v=0, lty = 2)


holein1[, 20] <- sapply(holein1[, 20], as.numeric)
holein1[, 22] <- sapply(holein1[, 22], as.numeric)
holein1[, 30] <- sapply(holein1[, 30], as.numeric)
mydata <- holein1[, c(2:6,20:30)]
chart.Correlation(mydata, histogram=TRUE, pch=19)


#normalisation
slope1<-as.data.frame(slope1[2,])
slope1
slope2<-as.data.frame(slope2[2,])
slope2
slope3<-as.data.frame(slope3[2,])
slope3
slope4<-as.data.frame(slope4[2,])
slope4
slope5<-as.data.frame(slope5[2,])
slope5

age
sexe
concentration
reads
RIN
cmon1<-as.data.frame(cmon1)
slope1
slope2
slope3
slope4
slope5

nowdata<-cbind(slope1,slope2,slope3,slope4,slope5,cmon1)#tableau avec intercepet et facteurs
attach(nowdata)
nowdata

sexe<-as.numeric(sexe)
sexe

#pour GB1,gene1 :

newdata<-`slope1[2, ]`[1]*age[1]+`slope2[2, ]`[1]*sexe[1]+`slope3[2, ]`[1]*reads[1]+`slope4[2, ]`[1]*concentration[1]+`slope5[2, ]`[1]*RIN[1]
newdata1<-`slope1[2, ]`[1]*age[2]+`slope2[2, ]`[1]*sexe[2]+`slope3[2, ]`[1]*reads[2]+`slope4[2, ]`[1]*concentration[2]+`slope5[2, ]`[1]*RIN[2]
newdata2<-`slope1[2, ]`[1]*age[3]+`slope2[2, ]`[1]*sexe[3]+`slope3[2, ]`[1]*reads[3]+`slope4[2, ]`[1]*concentration[3]+`slope5[2, ]`[1]*RIN[3]
newdata3<-`slope1[2, ]`[1]*age[4]+`slope2[2, ]`[1]*sexe[4]+`slope3[2, ]`[1]*reads[4]+`slope4[2, ]`[1]*concentration[4]+`slope5[2, ]`[1]*RIN[4]
newdata4<-`slope1[2, ]`[1]*age[5]+`slope2[2, ]`[1]*sexe[5]+`slope3[2, ]`[1]*reads[5]+`slope4[2, ]`[1]*concentration[5]+`slope5[2, ]`[1]*RIN[5]
newdata5<-`slope1[2, ]`[1]*age[6]+`slope2[2, ]`[1]*sexe[6]+`slope3[2, ]`[1]*reads[6]+`slope4[2, ]`[1]*concentration[6]+`slope5[2, ]`[1]*RIN[6]
newdata6<-`slope1[2, ]`[1]*age[7]+`slope2[2, ]`[1]*sexe[7]+`slope3[2, ]`[1]*reads[7]+`slope4[2, ]`[1]*concentration[7]+`slope5[2, ]`[1]*RIN[7]
newdata7<-`slope1[2, ]`[1]*age[8]+`slope2[2, ]`[1]*sexe[8]+`slope3[2, ]`[1]*reads[8]+`slope4[2, ]`[1]*concentration[8]+`slope5[2, ]`[1]*RIN[8]
newdata8<-`slope1[2, ]`[1]*age[9]+`slope2[2, ]`[1]*sexe[9]+`slope3[2, ]`[1]*reads[9]+`slope4[2, ]`[1]*concentration[9]+`slope5[2, ]`[1]*RIN[9]
newdata9<-`slope1[2, ]`[1]*age[10]+`slope2[2, ]`[1]*sexe[10]+`slope3[2, ]`[1]*reads[10]+`slope4[2, ]`[1]*concentration[10]+`slope5[2, ]`[1]*RIN[10]
newdata10<-`slope1[2, ]`[1]*age[11]+`slope2[2, ]`[1]*sexe[11]+`slope3[2, ]`[1]*reads[11]+`slope4[2, ]`[1]*concentration[11]+`slope5[2, ]`[1]*RIN[11]
newdata11<-`slope1[2, ]`[1]*age[12]+`slope2[2, ]`[1]*sexe[12]+`slope3[2, ]`[1]*reads[12]+`slope4[2, ]`[1]*concentration[12]+`slope5[2, ]`[1]*RIN[12]
newdata12<-`slope1[2, ]`[1]*age[13]+`slope2[2, ]`[1]*sexe[13]+`slope3[2, ]`[1]*reads[13]+`slope4[2, ]`[1]*concentration[13]+`slope5[2, ]`[1]*RIN[13]
newdata13<-`slope1[2, ]`[1]*age[14]+`slope2[2, ]`[1]*sexe[14]+`slope3[2, ]`[1]*reads[14]+`slope4[2, ]`[1]*concentration[14]+`slope5[2, ]`[1]*RIN[14]
newdata14<-`slope1[2, ]`[1]*age[15]+`slope2[2, ]`[1]*sexe[15]+`slope3[2, ]`[1]*reads[15]+`slope4[2, ]`[1]*concentration[15]+`slope5[2, ]`[1]*RIN[15]
newdata15<-`slope1[2, ]`[1]*age[16]+`slope2[2, ]`[1]*sexe[16]+`slope3[2, ]`[1]*reads[16]+`slope4[2, ]`[1]*concentration[16]+`slope5[2, ]`[1]*RIN[16]
newdata16<-cmon1[1]+`slope1[2, ]`[1]*age[17]+`slope2[2, ]`[1]*sexe[17]+`slope3[2, ]`[1]*reads[17]+`slope4[2, ]`[1]*concentration[17]+`slope5[2, ]`[1]*RIN[17]
newdata17<-cmon1[1]+`slope1[2, ]`[1]*age[18]+`slope2[2, ]`[1]*sexe[18]+`slope3[2, ]`[1]*reads[18]+`slope4[2, ]`[1]*concentration[18]+`slope5[2, ]`[1]*RIN[18]

error404<-cbind(newdata,newdata1,newdata2,newdata3,newdata4,newdata5,newdata6,newdata7,newdata8,newdata9,newdata10,newdata11,newdata12,newdata13,newdata14,newdata15,newdata16,newdata17)

error404  

names(error404)[1]<-"GB1"  
names(error404)[2]<-"GB2"
names(error404)[3]<-"GB3"  
names(error404)[4]<-"GB4"
names(error404)[5]<-"GB5"  
names(error404)[6]<-"GB6"
names(error404)[7]<-"GB7"  
names(error404)[8]<-"GB8"
names(error404)[9]<-"GB9"  
names(error404)[10]<-"GB10"
names(error404)[11]<-"GB11"  
names(error404)[12]<-"GB12"
names(error404)[13]<-"GB13"  
names(error404)[14]<-"GB14"
names(error404)[15]<-"GB15"  
names(error404)[16]<-"GB16"
names(error404)[17]<-"GB17"  
names(error404)[18]<-"GB18"

error405<-t(error404)
pca<-PCA(error405)








sapply(error404[],var)
range(sapply(error404[],var))
# maybe this range of variability is big in this context.
#Thus, we will use the correlation matrix
#For this, we must standardize our variables with scale() function:
iris.stand <- as.data.frame(scale(error404[,]))
sapply(iris.stand,sd) 

#If we use prcomp() function, we indicate 'scale=TRUE' to use correlation matrix
acp1 <- prcomp(iris.stand,scale=T)
#it is just the same that: prcomp(iris[,1:4],scale=T) and prcomp(iris.stand)
#similar with princomp(): princomp(iris.stand, cor=T)
acp1
summary(acp1)
#This gives us the standard deviation of each component, and the proportion of variance explained by each component.
#The stand
acp2<-as.data.frame(acp1$rotation)



#plot of variance of each PCA.
#It will be useful to decide how many principal components should be retained.
screeplot(acp1, type="lines",col=3)

acp1$rotation

summary(acp1)
names(summary(acp1))


pc<-head(pca)
head(acp2)

prcomp(acp1$rotation)
holein1<-merge(acp2, datavariable, by=0, all = T)


plot(acp1$rotation[,2], acp1$rotation[,3], pch = 19,  
     xlab="PC1 - 41.2%",ylab="PC2 - 18.4%")

abline(h=0, v=0, lty = 2)


holein1[, 18] <- sapply(holein1[, 18], as.numeric)
holein1[, 20] <- sapply(holein1[, 20], as.numeric)
holein1[, 28] <- sapply(holein1[, 28], as.numeric)
mydata <- holein1[, c(2:6,18:28)]
chart.Correlation(mydata, histogram=TRUE, pch=19)