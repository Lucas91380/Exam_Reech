############### EXERCICE 1 ###############

### Fonction de densité de l'énoncé
F_x <- function(x, a, b){
  return(ifelse(x >= b,(a*b^a)/(x^(a+1)),0))
}

#Plage x sur laquelle on va représenter nos fonctions
x <- seq(0,20,0.01)

#Représentation graphique de nos 3 fonction
plot(x, F_x(x=x,a=3,b=2), type='l', xlab='x', ylab='Densité', col='green')
lines(x, F_x(x=x,a=2,b=2), type='l', col='red')
lines(x, F_x(x=x,a=1,b=2), type = 'l', col='blue')
legend(15, 1.2, legend=c('a = 3', 'a = 2', 'a = 1'), col=c('green', 'red', 'blue'), lty=rep(1,3), cex=1)

#Fonction de simulation :
Simul_inv <- function(n, a, b){
  
  #Vecteur des u suivant une loi uniforme [0;1]
  u <- runif(n=n, min=0, max=1)
  
  #Fonction inverse de la fonction de répartition
  return((b/(-u+1)^(1/a)))
}

#Simulation d'un échantillon suivant notre densité
hist(Simul_inv(n=2000, a=1, b=2), xlim=c(0,20), breaks=5000, freq=FALSE)
lines(x, F_x(x=x,a=1,b=2), type = 'l', col='blue')

hist(Simul_inv(n=2000, a=2, b=2), xlim=c(0,20), breaks=500, freq=FALSE)
lines(x, F_x(x=x,a=2,b=2), type = 'l', col='red')

hist(Simul_inv(n=2000, a=3, b=2), xlim=c(0,20), breaks=200, freq=FALSE)
lines(x, F_x(x=x,a=3,b=2), type = 'l', col='green')



############### EXERCICE 2 ###############

#Fonction pour calculer notre Intégrale (Méthode uniforme)
Iunif <- function(x, m){
  
  #Tirage m valeurs suivants une loi uniforme (0;1)
  u = runif(m)  
  
  #Calcul de l'intégrale
  ihat = mean(1/sqrt(2*pi)*exp(-(-x/u)^2/2)*x/u^2) + 2*mean(1/sqrt(2*pi)*exp(-(u*x)^2/2)*x)
  
  #Calcul Vn
  vn = var(1/sqrt(2*pi)*exp(-(-x/u)^2/2)*x/u^2 + 2*1/sqrt(2*pi)*exp(-(u*x)^2/2)*x)
  
  return(list(ihat, vn))
}

#Fonction pour calculer notre Intégrale (Méthode normale)
Inorm <- function(x, m){
  
  #Tirage m valeurs suivants une loi normal (0;1)
  g = rnorm(m)
  
  #Calcul de l'intégrale
  ihat = mean(g<=x)
  
  #Calcul Vn
  vn = var(g<=x)
  
  return(list(ihat, vn))
}

#Paramètres de notre approximation
m <- 1000
x <- seq(0.1,3,0.1)

#Graphiques
plot(x, pnorm(x), type='l', ylim=c(0.5,1), ylab='F(x)')
lines(x, sapply(x,Iunif,m=1000)[1,], type='l', col='red')
lines(x, sapply(x,Inorm,m=1000)[1,], type='l', col='green')
legend(0,1, legend=c('Vrai F(x)', 'Estimation unif', 'Estimation norm'), col=c('black', 'red', 'green'), lty=rep(1,3), cex=1)

#Graphique des variances des estimateurs
plot(x, sapply(x,Iunif,m=1000)[2,], type='l', ylim=c(0,1) ,ylab='Var')
lines(x, sapply(x,Inorm,m=1000)[2,], type='l', col='red')
legend(0,1, legend=c('Uinf', 'Norm'), col=c('black', 'red'), lty=rep(1,2), cex=1)

m <- 1000

#Estimations
estim <- sapply(x,Iunif,m=m)

#Graphiques
plot(x, pnorm(x), type='l', ylim=c(0.5,1.2), ylab='F(x)', main = 'Estmimation par loi uniforme')
lines(x, unlist(estim[1,]), type='l', col='red')
lines(x, unlist(estim[1,])+1.96*sqrt(unlist(estim[2,])/m), type='l', col='red', lty=2)
lines(x, unlist(estim[1,])-1.96*sqrt(unlist(estim[2,])/m), type='l', col='red', lty=2)

m <- 1000

#Estimations
estim <- sapply(x,Inorm,m=m)

#Graphiques
plot(x, pnorm(x), type='l', ylim=c(0.5,1.2), ylab='F(x)', main = 'Estmimation par loi normale')
lines(x, unlist(estim[1,]), type='l', col='red')
lines(x, unlist(estim[1,])+1.96*sqrt(unlist(estim[2,])/m), type='l', col='red', lty=2)
lines(x, unlist(estim[1,])-1.96*sqrt(unlist(estim[2,])/m), type='l', col='red', lty=2)



############### EXERCICE 3 ###############

load("examen-dataset/spectra.Rdata")
library("FactoMineR")
library("factoextra")

set.seed(1)

## Vecteur de de longueur ncol(X) 
i <-  c(1:ncol(X))
## Sample aléatoire parmi les 546 spectres
j <- sample(nrow(X), 50, replace = FALSE)

## plot des spectres
plot(i,X[j[1],] , type="l", ylim = c(0.03, 0.09), xlab = "spectral channel",ylab="intensity", main = "Représentation de 50 spectres aléatoirements choisis")
for (k in j[2:50]) {
  lines(i,X[k,])
}

##ACP sur la matrice X 
res.pca <- PCA(X, graph = FALSE, scale.unit = TRUE)

##Extraction des valeurs propres
eig.val <- get_eigenvalue(res.pca)
var_exp <- eig.val[1:10,2]
var_cumul <- eig.val[2,3]


### Affichage des variances expliquées des 10 premiers axes
plot(var_exp, type = "b",xlab= "n-iéme composante",ylab="% de variance expliquée", main = "Proportion de variance expliquée par les 10 premières composantes")


### Plot des deux premiers axes sous formes de spectres (a finir)
plot(i,res.pca$var$coord[,1], type = "l",xlab = "spectral channel",ylab="intensity", main=" Spectres des deux premières composantes" );
lines(i,res.pca$var$coord[,2])

set.seed(13)
require(bootstrap)

n <- length(X[,1])
set.seed(1)
B <- 20
cor_boots_dim1 <- NULL
cor_boots_dim2 <- NULL
### ACP avec Bootstrap
var_exp_boots <- NULL
for(b in 1:B){
  ind = sample(c(1:n), n, replace = TRUE)
  acp_boots <- PCA(X[ind,], graph = FALSE, scale.unit = TRUE)
  eig.val_boots <- get_eigenvalue(acp_boots)
  var_exp_boots <- rbind(var_exp_boots,eig.val_boots[1:10,2])
  cor_boots_dim1 <- cbind(cor_boots_dim1,acp_boots$var$cor[,"Dim.1"])
  cor_boots_dim2 <- cbind(cor_boots_dim2,acp_boots$var$cor[,"Dim.2"])
}


# Calcul des intervalles de confiance
alpha = 0.05
Iperc <-NULL
b1 <- NULL; b2 <- NULL


# percentile method
for (j in 1:length(var_exp_boots[1,])) {
  q1 = quantile(var_exp_boots[,j], alpha/2)
  q2 = quantile(var_exp_boots[,j], 1-alpha/2)
  b2 <- rbind(b2,q2 )
  b1 <- rbind(b1, q1)
}
inter_conf <- cbind(b1, b2)


## Plot des la variance expliquée et de son intervalles de confiance
plot(inter_conf[,1], type = "l", col="red", main="Variance expliquée par composante", xlab="n-ième composante", ylab="% de variance expliquée", lty=2);lines(inter_conf[,2], col="red", type = "l", lty=2); lines(var_exp, type = "l")

hist(var_exp_boots[,"Dim.1"], col = "lightblue", border = "white",xlab="% de variance expliquée", main = "Histogramme de la variance expliquée \n par la 1ère composante estimée", ylab = "Fréquence")

abline(v = inter_conf[1,],lty = 2, lwd = 2, col = "red")

legend("topleft", c("red"), c("IC-perc"), col = c("red"),lty =2, lwd = 2, bg = "white")

abline(v = var_exp["Dim.1"], lwd = 2)

mtext("Valeur réelle", side = 1, at = var_exp["Dim.1"])

box()

hist(var_exp_boots[,"Dim.2"], col = "lightblue", border = "white",xlab="% de variance expliquée", main = "Histogramme de la variance expliquée \n par la 2ème composante estimée", ylab = "Fréquence")

abline(v = inter_conf[2,],lty = 2, lwd = 2, col = "red")

legend("topleft", c("red"), c("IC-perc"), col = c("red"),lty =2, lwd = 2, bg = "white")

abline(v = var_exp["Dim.2"], lwd = 2)

mtext("Valeur réelle", side = 1, at = var_exp["Dim.2"])

box()

##### Calcul des intervalles de confiance pour les cor des 2 premiers axes
alpha <- 0.05
Iperc <-NULL
b1 <- NULL; b2 <- NULL; b3 <- NULL; b4 <- NULL
inter_conf_cor <- NULL

##### percentile method
for (j in 1:length(cor_boots_dim2[,1])) {
  q1 <- quantile(cor_boots_dim1[j,], alpha/2)
  q2 <- quantile(cor_boots_dim1[j,], 1-alpha/2)
  q3 <- quantile(cor_boots_dim2[j,], alpha/2)
  q4 <- quantile(cor_boots_dim2[j,], 1-alpha/2)
  b1 <- rbind(b1,q1)
  b2 <- rbind(b2,q2 )
  b3 <- rbind(b3,q3 )
  b4 <- rbind(b4,q4 )
  
}
inter_conf_cor1 <- cbind(b1, b2)
inter_conf_cor2 <- cbind(b3, b4)


plot(i,res.pca$var$cor[,1], type = "l",xlab = "spectral channel",ylab="intensity", main=" Spectres des deux premières composantes avec interval de confiance" );
lines(i,res.pca$var$cor[,2])
lines(inter_conf_cor1[,1], col="red",lty =2)
lines(inter_conf_cor1[,2], col="red",lty =2)
lines(inter_conf_cor2[,1], col="blue",lty =2)
lines(inter_conf_cor2[,2], col="blue",lty =2)



############### EXERCICE 4 ###############

#Chargement des données
load("examen-dataset/permutation-dataset.Rdata")

#Affichage du graphique
plot(X, col = c(rep('black',100),rep('red', 100)),pch=16, xlab='x1', ylab='x2')

###Calcul de notre stat de Test

#Fonction qui trouve le voisin le plus proche
FindPpv <- function(ind){
  return(names(ind)[which(ind == min(ind, na.rm = TRUE))])
}

#Matrice des distances
distances <- as.matrix(dist(X))
diag(distances) <- NA

#Vecteur des voisins
ppv <- apply(distances,1,FindPpv)

#Calcul de la stat T
statT0 <- mean(c(as.numeric(ppv[1:100])<=100, as.numeric(ppv[101:200])>100))

#Paramètres b du nombre de répétition
b <- 1000
statT <- numeric(b)

#Boucle allant de 1 à b
for (i in 1:b){
  #Changement des positions des points au sein de X
  X_h0 <- X[sample(1:dim(X)[1], replace=F),]
  
  #Calcul de T
  distances <- as.matrix(dist(X_h0))
  diag(distances) <- NA
  ppv <- apply(distances,1,FindPpv)
  statT[i] <- mean(c(as.numeric(ppv[1:100])<=100, as.numeric(ppv[101:200])>100))
}

#Ajout de notre vrai T (observé sur les données) à la simulation
statT <- c(statT0,statT)

#Calcul de notre p-valeur
pval = mean(statT[-1] >= statT0)

#Affichage de l'histogramme de la distribution de nos T
hist(statT, xlab='T', main='Distribution de notre statistique T')
abline(v=statT0, col='red')

#Initialisation des paramètres
n <- seq(10,90,10)

#Fonction qui calcule fait tous les traitements de tirage etc et calcule la p-valeur directement
Pval_n <- function(n, t0, b){
  
  #Initialisation de statT
  statT <- numeric(b)
  
  #Boucle allant de 1 à b
  for (i in 1:b){
    
    #Tirage aléatoire de n valeurs de x dans chaque population
    subX <- rbind(X[sample(1:as.integer(dim(X)[1]/2), size=n, replace=F),], X[sample(as.integer(dim(X)[1]/2+1):dim(X)[1], size=n, replace=F),])
    
    #Changement des positions des points au sein de subX
    subX_h0 <- subX[sample(1:dim(subX)[1], replace=F),]
    
    #Calcul de T
    distances <- as.matrix(dist(subX_h0))
    diag(distances) <- NA
    ppv <- apply(distances,1,FindPpv)
    statT[i]<-mean(c(as.numeric(ppv[1:n])<=n, as.numeric(ppv[(n+1):(2*n)])>n))
  }
  statT <- c(t0,statT)
  return(mean(statT[-1] >= t0))
}

#Calcul de la p-valeur pour chaque n
pval_n <- sapply(n, Pval_n, t0=statT0, b=1000)

plot(x=n, y=pval_n, xlab='n', ylab='p-valeur', main='Evolution de la p-valeur en fonction de la taille de chaque groupe', type='b')

pval_n_10 <- sapply(n, Pval_n, t0=statT0, b=1000)
for(i in 1:9){
  pval_n_10 <- rbind(pval_n_10, sapply(n, Pval_n, t0=statT0, b=1000))
}

plot(x=n, y=pval_n_10[1,], xlab='n', ylab='p-valeur', main='Evolution de la p-valeur en fonction de la taille de chaque groupe', type='b', ylim=c(0,0.1))
for(i in 2:dim(pval_n_10)[1]){
  lines(x=n, y=pval_n_10[i,], xlab='n', ylab='p-valeur', main='Evolution de la p-valeur en fonction de la taille de chaque groupe', type='b', col=i) 
}
abline(h=0.05, col='red', lty=2)