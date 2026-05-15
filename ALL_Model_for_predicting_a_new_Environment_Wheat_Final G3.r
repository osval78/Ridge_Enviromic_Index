#################### Extending Finlay-Wilkinson regression with environmental covariates #########################################################################################################
#################### Initialization #########################################################################################################
#### Loading packages
#library(asreml)
rm(list = ls())
library(pls)
library(reshape2)  # Load the package
library(BGLR)
library(SKM)
library(reshape)
library(dplyr)
library(rrBLUP)

load("Wheat_G3_Data.RData")
ls()
EC=EC
Pheno=Pheno
Markers=Markers
GenoF=Geno
EC_YLD=EC_YLD
EC_TGW=EC_TGW
EC_GNO=EC_GNO
EC_HT=EC_HT
# Find the common lines present in all environments
common_lines <- Pheno %>%
  group_by(Line) %>%
  summarise(n_envs = n_distinct(Env)) %>%
  filter(n_envs == n_distinct(Pheno$Env)) %>%
  pull(Line)


# Filter the dataset to keep only the common lines
df_balanced <- Pheno %>%
  filter(Line %in% common_lines)

# Verify that all environments have the same lines
table(df_balanced$Env, df_balanced$Line)
PNames=colnames(df_balanced)
colnames(df_balanced)=c("Env","Genotype", PNames[-c(1:2)])
colnames(df_balanced)
final_geno_lines <- intersect(df_balanced$Genotype, rownames(GenoF)) %>% sort()
Geno <- GenoF[final_geno_lines, final_geno_lines]
dim(Geno)

lettuce0=merge(df_balanced,cbind(Env=EC[,1],scale(EC[,-1])),by="Env") 
lettuce0=lettuce0[order(lettuce0$Env,lettuce0$Genotype),]
head(lettuce0)
lettuceHT=merge(df_balanced,cbind(Env=EC_HT[,1],scale(EC_HT[,-1])),by="Env") 
lettuceHT=lettuceHT[order(lettuceHT$Env,lettuceHT$Genotype),]
head(lettuceHT)

lettuceGNO=merge(df_balanced,cbind(Env=EC_GNO[,1],scale(EC_GNO[,-1])),by="Env") 
lettuceGNO=lettuceGNO[order(lettuceGNO$Env,lettuceGNO$Genotype),]

lettuceTGW=merge(df_balanced,cbind(Env=EC_TGW[,1],scale(EC_TGW[,-1])),by="Env") 
lettuceTGW=lettuceTGW[order(lettuceTGW$Env,lettuceTGW$Genotype),]

lettuceYLD=merge(df_balanced,cbind(Env=EC_YLD[,1],scale(EC_YLD[,-1])),by="Env") 
lettuceYLD=lettuceYLD[order(lettuceYLD$Env,lettuceYLD$Genotype),]

# Convert relevant columns to numeric
# Convert relevant columns to numeric
numeric_vars <- c("DL", "GDD", "PTT", "PTR")


for (var in numeric_vars) {
  lettuce0[[var]] <- as.numeric(as.character(lettuce0[[var]]))
}
str(lettuce0)
CE=lettuce0[,-c(1:4)]
CE=as.matrix(CE)
K_CE=CE%*%t(CE)/ncol(CE)
L=t(chol(Geno))
lettuce0$Genotype=as.factor(lettuce0$Genotype)
X1 <- model.matrix(~0 + Genotype, data=lettuce0)  # Fixed effects
X1=X1%*%L

Z <- model.matrix(~Env - 1, data=lettuce0)  # Random environment effects
KE=Z%*%t(Z)/ncol(Z)
dim(Z)
KL=X1%*%t(X1)/ncol(X1)
KGE=KE*KL

ZE <- model.matrix(~0+Env, data=lettuce0)
colnames(ZE)
EC_S=scale(EC[,-1])
Phi=EC_S%*%t(EC_S)/ncol(EC_S)
KPhiE=ZE%*%Phi%*%t(ZE)
KGEP=KPhiE*KL
No_Interactions=5000
No_Burning=2000
############Leave one environment out (LOEO)
Traits_Names=colnames(lettuce0)[3:6]
Traits_Names
Summary_All=data.frame()
for (t in 1:length(Traits_Names)) {
#  t=1
  Trait=Traits_Names[t] 
  SKM::echo("*** Trait: %s ***", Trait)
  Traits_Names
  yy=lettuce0[,Trait]
  yy
  Summary=data.frame()
  Partitions=nrow(EC)
  for (i in 1:Partitions) {
#     i=3
    SKM::echo("\t*** Fold %s / %s ***", i, length(Partitions))
    set.seed(i)
    tst_set=i
    # environmental covariates matrices 
    Trn_data_Env=EC[-tst_set,-1]  
    Tst_data_Env=EC[tst_set,-1] 
    All_data_Env=EC[,-1]  
    Means_Trn=apply(Trn_data_Env,2,mean)
    SD_Trn=apply(Trn_data_Env,2,sd)
    C_Trn=scale(Trn_data_Env,center=Means_Trn,scale=SD_Trn) # C is identical to 'environment' data.frame but scaled
    C_Tst=scale(Tst_data_Env,center=Means_Trn,scale=SD_Trn) # C is identical to 'environment' data.frame but scaled
    C_All=scale(All_data_Env,center=Means_Trn,scale=SD_Trn)
    #################### PLS approach #########################################################################################################
    # Prepare the data
    Name_tst_Env=EC[tst_set,1] 
    pos_NA=which(lettuce0$Env==Name_tst_Env)
    rownames(lettuce0)=1:nrow(lettuce0)
    y <-yy  # Response variable
    yNA=y
    yNA[pos_NA]=NA
    ################Model M1 Conv without covariates#################
    ETAM1 <- list(
      FIXED1 = list(X =X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      RANDOM_ENV = list(K =KE, model = "RKHS"),  # Random ENV effect
      RANDOM_GE = list(K =KGE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM1 <- BGLR(y =yNA, ETA = ETAM1, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M1=mdl_bglrM1$yHat[pos_NA]
    Observed_tst=yy[pos_NA]
 
    ################Model M2 Conv with covariates#################
    ETAM2 <- list(
      FIXED1 = list(X =X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      K_CE=list(K =KGEP, model = "RKHS"),
      RANDOM_ENV = list(K =KPhiE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM2 <- BGLR(y =yNA , ETA = ETAM2, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M2=mdl_bglrM2$yHat[pos_NA]
    
    ############Model M3
    ETAM3 <- list(
      FIXED1 = list(X =X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      K_CE=list(K =KGEP, model = "RKHS"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM3 <- BGLR(y =yNA , ETA = ETAM3, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M3=mdl_bglrM3$yHat[pos_NA]
    
 ##########Models M4 and M5
    #lettuce0$FV <-mdl_bglrM1$yHat-XG%*%BetasG-mdl_bglr$mu
    lettuce0$FV <-mdl_bglrM2$ETA$K_CE$u
    
    # Convert FV into a matrix for SVD
    Y_matrix <- reshape::cast(data = lettuce0[, c("Genotype", "Env", "FV")], 
                              formula = Genotype ~ Env, value = "FV")[,-1]
    
    row.names(Y_matrix) <- levels(lettuce0$Genotype)
    #svd_result <-factanal(Y_matrix, factors=4)
    svd_result <- svd(Y_matrix)
    z1 <- svd_result$v[,1]  # First singular vector
    z2 <- svd_result$v[,2]  # Second singular vector

    lettuce=merge(x = lettuce0,y = data.frame(Env=unique(lettuce0$Env),z1=z1,z2=z2))
    lettuce=lettuce[order(lettuce$Env,lettuce$Genotype),]
    lettuce
    
    #############Model M2 with 1 latent covariate##
    X2 <- model.matrix(~0 + z1:X1, data=lettuce) 
    X3 <- model.matrix(~0 + z2:X1, data=lettuce) 

    # Define the priors for the random effects
    ETAM4 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
  #### Run the BGLR model
    mdl_bglrM4 <- BGLR(y = yNA, ETA = ETAM4, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M4=mdl_bglrM4$yHat[pos_NA]
  
    ###########Model M5
    # Define the priors for the random effects
    ETAM5 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      FIXED3 = list(X = X3, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM5 <- BGLR(y = yNA, ETA = ETAM5, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M5=mdl_bglrM5$yHat[pos_NA]
    
    ##############################################################################################
    #############Model M6 y M7 covariate## 
 
    lettuce0_Trn=lettuce0[-pos_NA,]
    lettuce0_Trn=droplevels(lettuce0_Trn)
    Y_matrix_Trn <- reshape::cast(data = lettuce0_Trn[, c("Genotype", "Env", Trait)], 
                                  formula = Env ~ Genotype, value = Trait)
    
    #################### PLS approach #########################################################################################################
    # run the PLS to get z1 and z2
    Y_data_trn=Y_matrix_Trn[,-1]
    Y_Trn=apply(Y_data_trn,1,mean,na.rm = TRUE)
    X_Trn <- as.matrix(C_Trn)
    X_Tst=as.matrix(C_Tst)
    Y_Trn=as.numeric(Y_Trn)
    
    pls_trn=mvr(Y_Trn~X_Trn, scale = FALSE, center = FALSE, ncomp = 3,method = "oscorespls") 
    
    Y_Hat1=c(predict(pls_trn,C_All,1))+mean(Y_Trn)
    Y_Hat2=c(predict(pls_trn,C_All,2))+mean(Y_Trn)
    
    Y_Hat1[-tst_set]=Y_Trn
    Y_Hat2[-tst_set]=Y_Trn
   
    z1=scale(Y_Hat1)
    z2=scale(Y_Hat2)
   
    lettuce=merge(x = lettuce0,y = data.frame(Env=unique(lettuce0$Env),z1=z1,z2=z2))
    lettuce=lettuce[order(lettuce$Env,lettuce$Genotype),]
    
    #############Model M2 with 1 latent covariate##
    X2 <- model.matrix(~0 + z1:X1, data=lettuce) 
    X3 <- model.matrix(~0 + z2:X1, data=lettuce) 
    
    # Define the priors for the random effects
    ETAM6 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM6 <- BGLR(y = yNA, ETA = ETAM6, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M6=mdl_bglrM6$yHat[pos_NA]
    ###########Model M7
    # Define the priors for the random effects
    ETAM7 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      FIXED3 = list(X = X3, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM7 <- BGLR(y = yNA, ETA = ETAM7, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M7=mdl_bglrM7$yHat[pos_NA]
    
    #########Models PLS M8 and M9######
    lettuce0_Trn=lettuce0[-pos_NA,]
    lettuce0_Trn=droplevels(lettuce0_Trn)
    Y_matrix_Trn <- reshape::cast(data = lettuce0_Trn[, c("Genotype", "Env", Trait)], 
                                  formula = Env ~ Genotype, value = Trait)
    
    #################### PLS Multivariate approach #########################################################################################################
    # run the PLS to get z1 and z2
    Y_data_trn=Y_matrix_Trn[,-1]
    Y_Trn=as.matrix(Y_data_trn)
    X_Trn=as.matrix(C_Trn)
    X_Tst=as.matrix(C_Tst)
    Y_Trn2=Y_Trn
    for (l in 1:ncol(Y_Trn)){
      # l=2
      Col_l=Y_Trn[,l]
      Mean_l=mean(Col_l, na.rm = TRUE)
      pos_CNA=which(is.na(Col_l))
      pos_CNA
      Y_Trn2[pos_CNA,l]=Mean_l
    }

    pls_trn=mvr(Y_Trn2~X_Trn,scale = FALSE,center =FALSE,ncomp =2,method = "oscorespls") 
    z1=pls_trn$scores[,1]
    z2=pls_trn$scores[,2]
    names(z1)=Y_matrix_Trn[,1]
    names(z2)=Y_matrix_Trn[,1]
    
    W_loading=pls_trn$loading.weights[,1:2]
    #####Matrix containing the regression coefficients
    P_loadings=pls_trn$loadings[,1:2]
    Pt_loadings=t(P_loadings)
    ####DtW inverse
    PtW_loadings_inv=solve(Pt_loadings%*%W_loading)
    ###R is equal to W multiplied by DtW_inv
    R=W_loading%*%PtW_loadings_inv
    ZS_Tst=X_Tst%*%R
    #### PLS1
    Z11=rep(0, Partitions)
    Z22=rep(0, Partitions)
    Z11[-tst_set]=z1
    Z11[tst_set]=ZS_Tst[1]
    Z22[-tst_set]=z2
    Z22[tst_set]=ZS_Tst[2]
    Z11=scale(Z11)
    Z22=scale(Z22)
    
    lettuce=merge(x = lettuce0,y = data.frame(Env=unique(lettuce0$Env),z1=Z11,z2=Z22))
    lettuce=lettuce[order(lettuce$Env,lettuce$Genotype),]
    
    #############Model M2 with 1 latent covariate##
    X2 <- model.matrix(~0 + z1:X1, data=lettuce) 
    X3 <- model.matrix(~0 + z2:X1, data=lettuce) 
    
    # Define the priors for the random effects
    ETAM8 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM8 <- BGLR(y = yNA, ETA = ETAM8, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M8=mdl_bglrM8$yHat[pos_NA]
    ###########Model M5
    # Define the priors for the random effects
    ETAM9 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      FIXED3 = list(X = X3, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM9 <- BGLR(y = yNA, ETA = ETAM9, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M9=mdl_bglrM9$yHat[pos_NA]
    
    ###########Model M10
    if(Trait=="YLD"){
    lettuce_Trait=lettuceYLD
    X2 <- model.matrix(~0 + PTR:X1, data=lettuce_Trait) 
    } else if (Trait=="GNO") {
      lettuce_Trait=lettuceGNO
      X2 <- model.matrix(~0 + PTR:X1, data=lettuce_Trait)
    }  else if (Trait=="TGW") {
      lettuce_Trait=lettuceTGW
      X2 <- model.matrix(~0 + PTT:X1, data=lettuce_Trait)
    } else {
      lettuce_Trait=lettuceHT
      X2 <- model.matrix(~0 + PTT:X1, data=lettuce_Trait)
    }
    
    # Define the priors for the random effects
    ETAM10 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM10 <- BGLR(y = yNA, ETA = ETAM10, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M10=mdl_bglrM10$yHat[pos_NA]
    
    # Define the priors for the random effects
    ETAM11 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR")
    )
    
    # Run the BGLR model
    mdl_bglrM11 <- BGLR(y = yNA, ETA = ETAM11, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M11=mdl_bglrM11$yHat[pos_NA]
        
##############################################################################################
    #############Model M12 covariate## 
    lettuce0_Trn=lettuce0[-pos_NA,]
    lettuce0_Trn=droplevels(lettuce0_Trn)
    Y_matrix_Trn <- reshape::cast(data = lettuce0_Trn[, c("Genotype", "Env", Trait)], 
                              formula = Env ~ Genotype, value = Trait)
    
    #################### Model Ridge Regression model M14 #########################################################################################################
    # run the PLS to get z1 and z2
    Y_data_trn=Y_matrix_Trn[,-1]
    Y_Trn=as.matrix(Y_data_trn)
    X_Trn=as.matrix(C_Trn)
    
    
    Y_Mean <- apply(Y_data_trn, 1, mean, na.rm = TRUE)
 
    SVD_X_All=svd(C_All)
    SVD_X_All
    UAll=SVD_X_All$u
    DAll=diag(SVD_X_All$d)
    
    X_NewAll=UAll%*%DAll
    X_NewAll
    
    pos_NAEnv=which(EC$Env==Name_tst_Env)
    pos_NAEnv
    X_New=X_NewAll[-pos_NAEnv,]
    # Fitting a GBLUP un-structured cov-matrices
    LP<-list(EC=list(X=X_New, model="BRR"))
    set.seed(123)

    fmUN<-BGLR(y=Y_Mean, ETA=LP, nIter=No_Interactions, burnIn=No_Burning,
               ,verbose=FALSE)
    
    #########Beta coeficients##################
    Mu=fmUN$mu
    BetasEC=fmUN$ETA$EC$b
    YY_Hat1=Mu+X_NewAll%*%BetasEC
    YY_Hat1
    
    YY_Hat1[-pos_NAEnv]=Y_Mean
   
    
    z1=scale(YY_Hat1[,1])
    z2=scale(YY_Hat1[,1])

    lettuce=merge(x = lettuce0,y = data.frame(Env=unique(lettuce0$Env),z1=z1,z2=z2,z3=z2))
    lettuce=lettuce[order(lettuce$Env,lettuce$Genotype),]
    
    X2 <- model.matrix(~0 + z1:X1, data=lettuce) 

    ##### Define the priors for the random effects
    ETAM12 <- list(
      FIXED1 = list(X = X1, model = "BRR"),  # Fixed effects with Bayesian Ridge Regression
      FIXED2 = list(X = X2, model = "BRR"),
      RANDOM_ENV = list(K =KE, model = "RKHS")  # Random ENV effect
    )
    
    # Run the BGLR model
    mdl_bglrM12 <- BGLR(y = yNA, ETA = ETAM12, nIter =No_Interactions, burnIn =No_Burning, verbose = FALSE)
    yhat_M12=mdl_bglrM12$yHat[pos_NA]

    #####Metricas M1 Convencional with covariates############
    MSE_M1=mse(yhat_M1,Observed_tst)
    COR_M1=cor(yhat_M1,Observed_tst,use="complete.obs")
    NRMSE_M1=nrmse(yhat_M1,Observed_tst)
    #####Metricas M2 GBLUP_RA1 with covariates############
    MSE_M2=mse(yhat_M2,Observed_tst)
    COR_M2=cor(yhat_M2,Observed_tst,use="complete.obs")
    NRMSE_M2=nrmse(yhat_M2,Observed_tst)
    #####Metricas M3 GBLUP_RA2 with covariates############
    MSE_M3=mse(yhat_M3,Observed_tst)
    COR_M3=cor(yhat_M3,Observed_tst,use="complete.obs")
    NRMSE_M3=nrmse(yhat_M3,Observed_tst)
    #####Metricas M4 ############
    MSE_M4=mse(yhat_M4,Observed_tst)
    COR_M4=cor(yhat_M4,Observed_tst,use="complete.obs")
    NRMSE_M4=nrmse(yhat_M4,Observed_tst)
    #####Metricas M5 ############
    MSE_M5=mse(yhat_M5,Observed_tst)
    COR_M5=cor(yhat_M5,Observed_tst,use="complete.obs")
    NRMSE_M5=nrmse(yhat_M5,Observed_tst) 
    #####Metricas M6 ############
    MSE_M6=mse(yhat_M6,Observed_tst)
    COR_M6=cor(yhat_M6,Observed_tst,use="complete.obs")
    NRMSE_M6=nrmse(yhat_M6,Observed_tst)
    #####Metricas M7 ############
    MSE_M7=mse(yhat_M7,Observed_tst)
    COR_M7=cor(yhat_M7,Observed_tst,use="complete.obs")
    NRMSE_M7=nrmse(yhat_M7,Observed_tst) 
    #####Metricas M8 ############
    MSE_M8=mse(yhat_M8,Observed_tst)
    COR_M8=cor(yhat_M8,Observed_tst,use="complete.obs")
    NRMSE_M8=nrmse(yhat_M8,Observed_tst) 
    #####Metricas M9 ############
    MSE_M9=mse(yhat_M9,Observed_tst)
    COR_M9=cor(yhat_M9,Observed_tst,use="complete.obs")
    NRMSE_M9=nrmse(yhat_M9,Observed_tst) 
    #####Metricas M10 ############
    MSE_M10=mse(yhat_M10,Observed_tst)
    COR_M10=cor(yhat_M10,Observed_tst,use="complete.obs")
    NRMSE_M10=nrmse(yhat_M10,Observed_tst) 
    #####Metricas M11 ############
    MSE_M11=mse(yhat_M11,Observed_tst)
    COR_M11=cor(yhat_M11,Observed_tst,use="complete.obs")
    NRMSE_M11=nrmse(yhat_M11,Observed_tst)
    #####Metricas M12 ############
    MSE_M12=mse(yhat_M12,Observed_tst)
    COR_M12=cor(yhat_M12,Observed_tst,use="complete.obs")
    NRMSE_M12=nrmse(yhat_M12,Observed_tst) 
    
    Summary_M1=data.frame(Method="M1",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M1,COR=COR_M1,NRMSE=NRMSE_M1,RE_MSE=MSE_M1/MSE_M12,RE_COR=COR_M12/COR_M1,RE_NRMSE=NRMSE_M1/NRMSE_M12)
    Summary_M2=data.frame(Method="M2",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M2,COR=COR_M2,NRMSE=NRMSE_M2,RE_MSE=MSE_M2/MSE_M12,RE_COR=COR_M12/COR_M2,RE_NRMSE=NRMSE_M2/NRMSE_M12)
    Summary_M3=data.frame(Method="M3",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M3,COR=COR_M3,NRMSE=NRMSE_M3,RE_MSE=MSE_M3/MSE_M12,RE_COR=COR_M12/COR_M3,RE_NRMSE=NRMSE_M3/NRMSE_M12)
    Summary_M4=data.frame(Method="M4",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M4,COR=COR_M4,NRMSE=NRMSE_M4,RE_MSE=MSE_M4/MSE_M12,RE_COR=COR_M12/COR_M4,RE_NRMSE=NRMSE_M4/NRMSE_M12) 
    Summary_M5=data.frame(Method="M5",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M5,COR=COR_M5,NRMSE=NRMSE_M5,RE_MSE=MSE_M5/MSE_M12,RE_COR=COR_M12/COR_M5,RE_NRMSE=NRMSE_M5/NRMSE_M12)
    Summary_M6=data.frame(Method="M6",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M6,COR=COR_M6,NRMSE=NRMSE_M6,RE_MSE=MSE_M6/MSE_M12,RE_COR=COR_M12/COR_M6,RE_NRMSE=NRMSE_M6/NRMSE_M12)
    Summary_M7=data.frame(Method="M7",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M7,COR=COR_M7,NRMSE=NRMSE_M7,RE_MSE=MSE_M7/MSE_M12,RE_COR=COR_M12/COR_M7,RE_NRMSE=NRMSE_M7/NRMSE_M12)
    Summary_M8=data.frame(Method="M8",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M8,COR=COR_M8,NRMSE=NRMSE_M8,RE_MSE=MSE_M8/MSE_M12,RE_COR=COR_M12/COR_M8,RE_NRMSE=NRMSE_M8/NRMSE_M12)
    Summary_M9=data.frame(Method="M9",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M9,COR=COR_M9,NRMSE=NRMSE_M9,RE_MSE=MSE_M9/MSE_M12,RE_COR=COR_M12/COR_M9,RE_NRMSE=NRMSE_M9/NRMSE_M12)
    Summary_M10=data.frame(Method="M10",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M10,COR=COR_M10,NRMSE=NRMSE_M10,RE_MSE=MSE_M10/MSE_M12,RE_COR=COR_M12/COR_M10,RE_NRMSE=NRMSE_M10/NRMSE_M12)
    Summary_M11=data.frame(Method="M11",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M11,COR=COR_M11,NRMSE=NRMSE_M11,RE_MSE=MSE_M11/MSE_M12,RE_COR=COR_M12/COR_M11,RE_NRMSE=NRMSE_M11/NRMSE_M12)
    Summary_M12=data.frame(Method="M12",Trait=Trait,Env=Name_tst_Env,MSE=MSE_M12,COR=COR_M12,NRMSE=NRMSE_M12,RE_MSE=MSE_M12/MSE_M12,RE_COR=COR_M12/COR_M12,RE_NRMSE=NRMSE_M12/NRMSE_M12)

    Summary_i=rbind(Summary_M1,Summary_M2,Summary_M3,Summary_M4,Summary_M5,Summary_M6,Summary_M7,Summary_M8,Summary_M9,Summary_M10,Summary_M11,Summary_M12)
    Summary=rbind(Summary,Summary_i)
    Summary
  }
  Summary_All=rbind(Summary_All,Summary)
}
Summary_All
write.csv(Summary_All,file="All_Methods_Wheat_Final_G3_Final_Nuevo.csv")



