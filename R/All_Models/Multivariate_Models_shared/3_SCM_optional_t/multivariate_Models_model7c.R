################################################################################
################################################################################
######               Fitting Multivariate Models. Model 7c                ######
################################################################################
################################################################################
library(INLA)
inla.setOption(inla.mode = "compact")
library(spdep)
library(maptools)
library(MASS)

########################################
## Data organization for INLA         ##
########################################
## Load data and SpatialPolygonsDataFrame
carto <- st_read("../../../Data/Carto/carto_gb.shp")
carto <- carto[order(carto$Code), ]

## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
carto_use <- sf::st_as_sf(carto)

# adjacency matrix
adj_gb <- as.matrix(read.table("../../../Data/Carto/adj_gb.txt"))

carto.nb <- mat2listw(adj_gb)$neighbours
W <- as(nb2mat(carto.nb, style = "B"), "Matrix")

spdep::nb2INLA("carto_nb.graph", carto.nb)
g <- INLA::inla.read.graph("carto_nb.graph")


## Load data 
load("../../../Data/Data_Pancreas_male.rda")
data <- Data.Pancreas

n <- length(unique(data$Code))
t <- length(unique(data$ID_time))
J <- length(unique(data$health_outcome))


###Index for INLA. Important Note: The order of the data is crucial due to the 
# shared spatial effect. It is necessary to designate areas from 1 to n for 
# incidence and areas from (n+1) to (2*n) for mortality.
data$ID_type<-rep(c(1,2),each=n*t)

data$alpha1 <- rep(1:0,each=n*t)
data$alpha2 <- rep(0:1,each=n*t)

data$ID_area <- rep(1:n,t*2)
data$ID_area1 <- data$ID_area+n*(data$ID_time-1) + c(rep(0,n*t),rep(n*t,n*t))

data$ID_area <- c(rep(1:n,t),rep((n+1):(2*n),t))
data$ID_unst2 <- c(rep(1:n,t),rep(1:n,t))

data$ID_time1 <- c(data$ID_time[which(data$ID_type==1)],rep(NA,n*t))
data$ID_time2 <- c(rep(NA,n*t),data$ID_time[which(data$ID_type==2)])


########################################
##  Constraints                ##
########################################
##Q
Q_xi <- matrix(0, g$n, g$n)
for (i in 1:g$n){
  Q_xi[i,i]=g$nnbs[[i]]
  Q_xi[i,g$nbs[[i]]]=-1
}

##RW1
D1 <- diff(diag(t), differences=1)
Q_gammaRW1 <- t(D1)%*%D1

## constraints
R_1_2 <- kronecker(Q_gammaRW1,diag(n))
r_def_1_2 <- n
A_constr_1_2<- kronecker(matrix(1,1,t), diag(n))

R_1_3 <- kronecker(diag(t),Q_xi)
r_def_1_3 <- t
A_constr_1_3<- kronecker(diag(t),matrix(1,1,n))

R_1_4 <- kronecker(Q_gammaRW1,Q_xi)
r_def_1_4 <- n+t-1
A.1.1 <- kronecker(matrix(1,1,t),diag(n))
A.1.2 <- kronecker(diag(t),matrix(1,1,n))
A_constr_1_4 <- rbind(A.1.1,A.1.2)


########################################
##  Priors                ##
########################################
#### UNIFORM PRIORS
sdunif="expression:

logdens=-log_precision/2;

return(logdens)"



################################################################################
######                            Model Fitting                           ######
################################################################################
########################################
##      Type I           ##
########################################
source("./inla_rgeneric_scm_change_typeI.R")

s = 3 #Number of different scaling parameters (l in the paper)
N = c(3, 3, 3) #Number of times each scaling parameter is repeated (m_l in the paper)
SCM.model.typeI <- inla.rgeneric.define(inla.rgeneric.SCM, debug = FALSE, 
                                        k = t, s = s, N = N,
                                        W = W, initial.values= c(4,rep(0,s)))


A.constr.s4 <- kronecker(diag(2), matrix(1,1,n*t))

formula = counts~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
    scale.model = T, constr = TRUE) +
  f(ID_unst2, model="iid", hyper=list(prec=list(prior=sdunif)), replicate = ID_type) +
  f(ID_time1, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_time2, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_area1, model= SCM.model.typeI, constr = FALSE, extraconstr=list(A=A.constr.s4, e=rep(0,2)))

shared_rw1_indep_t1 =  inla(formula,
                            family = "poisson",
                            data = data,
                            E=Population,
                            control.compute=list(dic=TRUE,
                                                 cpo=TRUE,
                                                 waic=TRUE,
                                                 hyperpar=TRUE,
                                                 config = TRUE,
                                                 return.marginals.predictor=TRUE),
                            control.inla=list(strategy="simplified.laplace",
                                              verbose = T,
                                              numint.maxfeval= 100000),
                            control.predictor = list(link=1,compute = TRUE))

########################################
##      Type II           ##
########################################
source("./inla_rgeneric_scm_change_typeII.R")

s = 3
N = c(3, 3, 3)
SCM.model <- inla.rgeneric.define(inla.rgeneric.SCM, debug = FALSE, 
                                  k = t, s = s, N = N,
                                  W = W, initial.values= c(4,rep(0,s)))


A.constr.s4 <- kronecker(diag(2),kronecker(matrix(1,1,t),diag(n)))#Aldatu beharra??? Efecktu espaziala eta tenporala sartzean

formula = counts~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
    scale.model = T, constr = TRUE) +
  f(ID_unst2, model="iid", hyper=list(prec=list(prior=sdunif)), replicate = ID_type) +
  f(ID_time1, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_time2, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_area1, model= SCM.model, constr = FALSE, extraconstr=list(A=A.constr.s4, e=rep(0,2*(n))))

shared_rw1_indep_t2 =  inla(formula,
                            family = "poisson",
                            data = data,
                            E=Population,
                            control.compute=list(dic=TRUE,
                                                 cpo=TRUE,
                                                 waic=TRUE,
                                                 hyperpar=TRUE,
                                                 config = TRUE,
                                                 return.marginals.predictor=TRUE),
                            control.inla=list(strategy="simplified.laplace",
                                              verbose = T,
                                              numint.maxfeval= 100000),
                            control.predictor = list(link=1,compute = TRUE))


########################################
##      Type III           ##
########################################
source("./inla_rgeneric_scm_change_typeIII.R")

s = 3
N = c(3, 3, 3)
SCM.model <- inla.rgeneric.define(inla.rgeneric.SCM, debug = FALSE, 
                                  k = t, s = s, N = N,
                                  W = W, initial.values= c(4,rep(0,s)))


A.constr.s4 <- kronecker(diag(2*t), matrix(1,1,n))

formula = counts~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
    scale.model = T, constr = TRUE) +
  f(ID_unst2, model="iid", hyper=list(prec=list(prior=sdunif)), replicate = ID_type) +
  f(ID_time1, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_time2, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_area1, model= SCM.model, constr = FALSE, extraconstr=list(A=A.constr.s4, e=rep(0,2*t)))

shared_rw1_indep_t3 =  inla(formula,
                            family = "poisson",
                            data = data,
                            E=Population,
                            control.compute=list(dic=TRUE,
                                                 cpo=TRUE,
                                                 waic=TRUE,
                                                 hyperpar=TRUE,
                                                 config = TRUE,
                                                 return.marginals.predictor=TRUE),
                            control.inla=list(strategy="simplified.laplace",
                                              verbose = T,
                                              numint.maxfeval= 100000),
                            control.predictor = list(link=1,compute = TRUE))



########################################
##      Type IV           ##
########################################
source("./inla_rgeneric_scm_change_typeIV.R")

s = 3
N = c(3, 3, 3)
SCM.model <- inla.rgeneric.define(inla.rgeneric.SCM, debug = FALSE, 
                                  k = t, s = s, N = N,
                                  W = W, initial.values= c(4,rep(0,s)))


A.constr.s1 <- kronecker(matrix(1,1,t),diag(n))#Aldatu beharra??? Efecktu espaziala eta tenporala sartzean
A.constr.s2 <- kronecker(diag(t), matrix(1,1,n))  #Aldatu beharra???
A.constr <- rbind(A.constr.s1,A.constr.s2)
A.constr.s4 <- kronecker(diag(2), A.constr)

formula = counts~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g, hyper=list(prec=list(prior=sdunif)),
    scale.model = T, constr = TRUE) +
  f(ID_unst2, model="iid", hyper=list(prec=list(prior=sdunif)), replicate = ID_type) +
  f(ID_time1, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_time2, model="rw1", hyper=list(prec=list(prior=sdunif)), constr = TRUE) +
  f(ID_area1, model= SCM.model, constr = FALSE, extraconstr=list(A=A.constr.s4, e=rep(0,2*(n+t))))


shared_rw1_indep_t4 = tryCatch({ inla(formula,
                                      family = "poisson",
                                      data = data,
                                      E=Population,
                                      control.compute=list(dic=TRUE,
                                                           cpo=TRUE,
                                                           waic=TRUE,
                                                           hyperpar=TRUE,
                                                           config = TRUE,
                                                           return.marginals.predictor=TRUE),
                                      control.inla=list(strategy="simplified.laplace",
                                                        verbose = T,
                                                        numint.maxfeval= 100000),
                                      control.predictor = list(link=1,compute = TRUE))},
                               error = function(msg) {
                                 print(msg)
                                 return(NULL)
                               })




################################################################################
######                      Model selection criteria                      ######
################################################################################
###Table for results:
tab <- matrix(NA, nrow=4, ncol=3)
index1 <- c("indep")
index2 <- c("t1","t2","t3","t4")
for (i in 1:length(index1)) {
  for (l in 1:length(index2)) {
    if (eval(parse(text = paste0("length(shared_rw1_",index1[i],"_",index2[l],")!=0")))){
      eval(parse(text = paste0("tab[(i-1)*length(index2)+l,] <- c(shared_rw1_",index1[i],"_",index2[l],"$dic$dic,
                           shared_rw1_",index1[i],"_",index2[l],"$waic$waic
                           ,sum(-log(shared_rw1_",index1[i],"_",index2[l],"$cpo$cpo)))")))
      
    }
    else{
      tab[(i-1)*length(index2)+l,] <- c(NA,NA,NA)
    }
  }
}
tab <- cbind(rep(index1,each=length(index2)),rep(index2,length(index1)),tab)
tab <- rbind(c(" "," ","DIC","WAIC","LS"), tab)

latex_table <- xtable::xtable(tab, caption='Model selection criteria. Model 7c', digits=3)
xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))


