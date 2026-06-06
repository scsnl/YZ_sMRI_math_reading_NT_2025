setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
library(CCA)
library(CCP) # for cca statistical test
library(permute)
library(readxl)
library(R.matlab) 
library(psych) # for pca
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggseg)

###########################################
### Identify extreme outliers using IQR ###
###########################################

identify_outliers <- function(column) {
  Q1 <- quantile(column, 0.25)
  Q3 <- quantile(column, 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 3 * IQR
  upper_bound <- Q3 + 3 * IQR
  which(column < lower_bound | column > upper_bound)
}


#################
### load data ###
#################

load_data <- function(bn_atlas_file, beh_file, brain_file){
  # roi labels
  bn_atlas = read_excel(bn_atlas_file, sheet=1)
  roi.names = bn_atlas$Description[1:218]
  
  # behavior
  beh = read.csv(beh_file)
  
  # brain
  brain = read.csv(brain_file)
  pids = brain[,1]
  brain = as.matrix(brain[,-c(1:3)])
  colnames(brain) = roi.names
  numsub = dim(brain)[1]
  numroi = dim(brain)[2]
  
  # output
  out = list(bn_atlas=bn_atlas, roi.names=colnames(brain),
             beh=beh, brain=brain, pids = pids, numsub=numsub, numroi=numroi)
  
  out
}


###########
### PCA ###
###########
mypca <- function(data, expvar){
  my.pca = prcomp(data, center=T, scale=F)
  
  var_explained = my.pca$sdev^2/sum(my.pca$sdev^2)
  var_cum = cumsum(var_explained)
  numPC = which(var_cum > expvar)[1]
  
  out = list(center = my.pca$center, scale = my.pca$scale, sdev = my.pca$sdev,
             rotation = my.pca$rotation, x = my.pca$x, 
             numPC = numPC)
  
  out
}


##########################
##### CCA function ######
##########################

mycca_noscale <- function(X,Y,nperm){
  
  epsilon = sqrt(.Machine$double.eps)
  
  X = scale(X, center = TRUE, scale = FALSE)
  Y = scale(Y, center = TRUE, scale = FALSE)
  x_mtx = X
  y_mtx = Y
  
  Y.c <- scale(Y, center = TRUE, scale = FALSE)
  X.c <- scale(X, center = TRUE, scale = FALSE)
  S11 <- cov(Y)
  S22 <- cov(X)
  S12 <- cov(Y,X)
  S11.chol <- chol(S11)
  S11.chol.inv <- solve(S11.chol)
  S22.chol <- chol(S22)
  S22.chol.inv <- solve(S22.chol)
  # K summarizes the correlation structure between the two sets of variables
  K <- t(S11.chol.inv) %*% S12 %*% S22.chol.inv
  K.svd <- svd(K)
  Eigenvalues <- K.svd$d^2
  axenames <- paste("CanAxis",seq_along(K.svd$d),sep="")
  U <- K.svd$u
  V <- K.svd$v
  A <- S11.chol.inv %*% U # raw canonical coefficients
  B <- S22.chol.inv %*% V
  Cy <- (Y %*% A)
  Cx <- (X %*% B)
  ## Compute the 'Biplot scores of Y and X variables' a posteriori
  corr.Y.Cy <- cor(Y.c, Cy)  # To plot Y in biplot in space Y
  corr.Y.Cx <- cor(Y.c, Cx)  # Available for plotting Y in space of X
  corr.X.Cy <- cor(X.c, Cy)  # Available for plotting X in space of Y
  corr.X.Cx <- cor(X.c, Cx)  # To plot X in biplot in space X
  ## Add names
  colnames(corr.Y.Cy) <- colnames(corr.Y.Cx) <- axenames
  colnames(corr.X.Cy) <- colnames(corr.X.Cx) <- axenames
  colnames(A) <- colnames(B) <- axenames
  rownames(A) <- colnames(Y)
  rownames(B) <- colnames(X)
  # standardized canonical coefficient
  Astd <- diag(sqrt(diag(cov(Y)))) %*% A
  Bstd <- diag(sqrt(diag(cov(X)))) %*% B
  
  # Compute Pillai's trace = sum of the canonical eigenvalues
  #                        = sum of the squared canonical correlations
  S11.inv <- S11.chol.inv %*% t(S11.chol.inv)
  S22.inv <- S22.chol.inv %*% t(S22.chol.inv)
  gross.mat <- S12 %*% S22.inv %*% t(S12) %*% S11.inv
  PillaiTrace <- sum(diag(gross.mat))
  
  n = nrow(X)
  pp = ncol(y_mtx)
  qq = ncol(x_mtx)
  s = min(pp,qq)
  df1 = max(pp,qq)
  df2 = (n - max(pp,qq) - 1)
  
  Fval  <- (PillaiTrace*df2)/((s-PillaiTrace)*df1)
  Fref <- Fval
  p.Pillai <- pf(Fval, s*df1, s*df2, lower.tail=FALSE)
  
  # permutation tests 
  nperm = nperm
  perm.idx = shuffleSet(n,nperm) # a row = a permutation including n (number of observations) indices
  
  # permutation test for pillai (model fit)
  Pillai.perm <- numeric(nperm)
  F.Pillai.perm <- numeric(nperm)

  for(perm in c(1:nperm)){
    idx = perm.idx[perm,]
    S12.per <- cov(Y[idx,], X)
    gross.mat.per <- S12.per %*% S22.inv %*% t(S12.per) %*% S11.inv
    Pillai.per <- sum(diag(gross.mat.per))
    Fper  <- (Pillai.per*df2)/((s-Pillai.per)*df1)
    
    Pillai.perm[perm] <- Pillai.per
    F.Pillai.perm[perm] <- Fper  
  }
  pillai.p.perm <- (sum(F.Pillai.perm >= (Fref - epsilon)) + 1) / (nperm + 1)
  
  
  # permutation test for canonical correlation (with CCA rerun)
  r.actual = diag(cor(Cx,Cy))
  
  cancor.p.perm = c()
  df.rrand = c()
  for(perm in c(1:nperm)){
    idx = perm.idx[perm,]
    tmp = cc(X,Y[idx,])
    r.rand = tmp$cor
    df.rrand = rbind(df.rrand,r.rand)
  }
  
  df.ractual = do.call("rbind", replicate(nperm, t(r.actual), simplify = FALSE))
  # r.perm = df.rrand > df.ractual
  # cancor.p.perm = colSums(r.perm) / dim(r.perm)[1]
  # cancor.p.perm.adjust = p.adjust(cancor.p.perm,method="fdr",n=length(r.actual))
  
  r.perm = df.rrand >= df.ractual
  cancor.p.perm = (colSums(r.perm) + 1) / (dim(r.perm)[1] + 1)
  cancor.p.perm.adjust = p.adjust(cancor.p.perm, method = "fdr", n = length(r.actual))
  
  # output
  out = list(Pillai=PillaiTrace, Eigenvalues=Eigenvalues, cancor=K.svd$d,
             xcoef.raw = B, ycoef.raw = A,
             xcoef.std = Bstd, ycoef.std = Astd,
             pillai.p = p.Pillai, pillai.p.perm = pillai.p.perm, 
             Pillai.perm = Pillai.perm, F.Pillai.perm = F.Pillai.perm,
             Cy=Cy, Cx=Cx,
             corr.Y.Cy=corr.Y.Cy, corr.X.Cx=corr.X.Cx, 
             corr.Y.Cx=corr.Y.Cx, corr.X.Cy=corr.X.Cy,
             cancor.rand=df.rrand, 
             cancor.p.perm=cancor.p.perm, cancor.p.perm.adjust=cancor.p.perm.adjust,
             nperm=nperm, perm.idx=perm.idx)
  
  out
}


################
### Plotting ###
################
myplot_YCy <- function(res, output_path, prefix, condition, postfix){
  ## Correlation between behavioral variables (Y) and modes (i.e. loadings) ##
  Y.Cy.ld = melt(res$corr.Y.Cy)
  df = Y.Cy.ld 
  yy = max(abs(min(df$value)),abs(max(df$value)))

  # heatmap color & value overlay
  ggplot() +
    geom_tile(df, mapping = aes(Var2, Var1, fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1*yy,yy)) +
    geom_text(df,mapping = aes(Var2, Var1,
                               label = round(value, digit = 3)),color="black") +
    labs(x="Modes", y="Behaviors") +
    theme_classic() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                            axis.text.y=element_text(size=9, angle=0, vjust=0.3),
                            plot.title=element_text(size=11))
  plt.name = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_YCy_heatmap.eps',sep="")
  ggsave(plt.name,units="in", width=4, height=5)
}

##
myplot_YcoefStd <- function(res, output_path, prefix, condition, postfix){
  ## Standardized weights for behavioral variables (Y) ## 
  rownames(res$ycoef.std) = rownames(res$corr.Y.Cy)
  ycoef.std = melt(res$ycoef.std)
  df = ycoef.std 
  yy = max(abs(min(df$value)),abs(max(df$value)))
  
  # heatmap color & value overlay
  ggplot() +
    geom_tile(df, mapping = aes(Var2, Var1, fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1*yy,yy)) +
    geom_text(df,mapping = aes(Var2, Var1,
                               label = round(value, digit = 3)),color="black") +
    labs(x="Modes", y="Behaviors") +
    theme_classic() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                            axis.text.y=element_text(size=9, angle=0, vjust=0.3),
                            plot.title=element_text(size=11))
  plt.name = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_Ystdcoef_heatmap.eps',sep="")
  ggsave(plt.name,units="in", width=4, height=5)
}


myplot_XCx <- function(res, output_path, prefix, condition, postfix){
  ## Correlation between behavioral variables (Y) and modes (i.e. loadings) ##
  X.Cx.ld = melt(res$corr.X.Cx)
  df = X.Cx.ld 
  yy = max(abs(min(df$value)),abs(max(df$value)))
  
  # heatmap color & value overlay
  ggplot() +
    geom_tile(df, mapping = aes(Var2, Var1, fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1*yy,yy)) +
    geom_text(df,mapping = aes(Var2, Var1,
                               label = round(value, digit = 3)),color="black") +
    labs(x="Modes", y="Brain") +
    theme_classic() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                            axis.text.y=element_text(size=9, angle=0, vjust=0.3),
                            plot.title=element_text(size=11))
  plt.name = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_XCx_heatmap.eps',sep="")
  ggsave(plt.name,units="in", width=4, height=5)
}

##
myplot_XcoefStd <- function(res, output_path, prefix, condition, postfix){
  ## Standardized weights for behavioral variables (Y) ## 
  rownames(res$xcoef.std) = rownames(res$corr.X.Cx)
  xcoef.std = melt(res$xcoef.std)
  df = xcoef.std 
  yy = max(abs(min(df$value)),abs(max(df$value)))
  
  # heatmap color & value overlay
  ggplot() +
    geom_tile(df, mapping = aes(Var2, Var1, fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1*yy,yy)) +
    geom_text(df,mapping = aes(Var2, Var1,
                               label = round(value, digit = 3)),color="black") +
    labs(x="Modes", y="Brain") +
    theme_classic() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                            axis.text.y=element_text(size=9, angle=0, vjust=0.3),
                            plot.title=element_text(size=11))
  plt.name = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_Xstdcoef_heatmap.eps',sep="")
  ggsave(plt.name,units="in", width=4, height=5)
}


##
myplot_cca_scatter <- function(res, output_path, prefix, condition, postfix){
  ## plot scatter plot (canonical scores) ##

  rownames(res$ycoef.std) = rownames(res$corr.Y.Cy)
  ycoef.std = melt(res$ycoef.std)
  df = ycoef.std 
  type = "colorBystdcoef"
  
  for(i in 1:dim(res$Cy)[2]){
    col_id = which(df$Var2 == paste("CanAxis",i, sep=""))
    idx = which(abs(df$value[col_id]) == max(abs(df$value[col_id])))
    
    print(df$Var1[idx])
    
    Cy = res$Cy[,i] # behavior
    Cx = res$Cx[,i] # brain
    df_plot = data.frame(Cx = Cx, Cy=Cy, Cog=Y[,idx])
    
    ggplot(df_plot, aes(x=Cx, y=Cy, col=Cog)) +
      geom_point() + 
      scale_colour_gradientn(colours = c("blue", "cyan", "yellow", "red"),
                             name = df$Var1[idx]) +
      labs(x=paste("Mode",i,"(Brain)",sep=" "), y=paste("Mode",i,"(Behavior)",sep=" ")) +
      theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=12, angle=0, vjust=0.3),
                         plot.title=element_text(size=12))
    plt.name = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_CxCy_scatter_mode', i, '_', type, '.eps',sep="")
    ggsave(plt.name,units="in", width=5, height=4)
  }
}


####################################################################
## write CCA brain loadings and std coef to *.mat for brain plots ##
####################################################################
write_weights <- function(res, Xorig, Y, apply_pca, output_path, prefix, condition, postfix){
  if(apply_pca == 1){
    # std coefs mapping from PCs to original ct rois
    roi.xcoef.std.allmodes = c()
    roi.xcoef.raw.allmodes = c()
    
    for(i in 1:dim(res$xcoef.std)[2]){
      roi.xcoef.std = t(res$xcoef.std[,i]) %*% t(my.pca$rotation[,c(1:my.pca$numPC)])
      roi.xcoef.raw = t(res$xcoef.raw[,i]) %*% t(my.pca$rotation[,c(1:my.pca$numPC)])
      if(i==1){
        roi.xcoef.std.allmodes = t(roi.xcoef.std)
        roi.xcoef.raw.allmodes = t(roi.xcoef.raw)
      }else{
        roi.xcoef.std.allmodes = cbind(roi.xcoef.std.allmodes,t(roi.xcoef.std))
        roi.xcoef.raw.allmodes = cbind(roi.xcoef.raw.allmodes,t(roi.xcoef.raw))
      }
    }
    roi.xcoef.std.allmodes = as.data.frame(roi.xcoef.std.allmodes)
    roi.xcoef.raw.allmodes = as.data.frame(roi.xcoef.raw.allmodes)
    
    # yz added - compute loadings (correlation between each mode and original variables)
    BrainCx = cor(Xorig, res$Cx)
    BehCy = cor(Y, res$Cy)
    
    # write std coef and raw coef and loadings
    file.name = paste(output_path,'/', prefix,'_',condition, '_',postfix,'_coef.mat',sep="")
    writeMat(file.name, roi.xcoef.std = roi.xcoef.std.allmodes, roi.xcoef.raw = roi.xcoef.raw.allmodes, 
             ycoef.std = res$ycoef.std, ycoef.raw = res$ycoef.raw, 
             BrainCx = BrainCx, BehCy = BehCy,
             yrownames=rownames(res$corr.Y.Cy))
    
    
    df_all = cbind(roi.xcoef.std.allmodes,BrainCx)
    
    num_modes <- ncol(df_all) / 2
    # Generate names
    xcoef_names <- paste0("xcoef_std_m", 1:num_modes)
    xloading_names <- paste0("xloading_m", 1:num_modes)
    # Assign column names
    colnames(df_all) <- c(xcoef_names, xloading_names)
    
    file.name = paste(output_path,'/', prefix,'_',condition, '_',postfix,'_coef.csv',sep="")
    write.csv(df_all, file.name)
    
  }else{
    ## Correlation between brain variables (original ct measures rather than PCs) 
    ## and modes (i.e. loadings) ##
    
    file.name = paste(output_path,'/', prefix,'_',condition, '_',postfix,'_corrXCx.mat',sep="")
    writeMat(file.name, BrainCx = res$corr.X.Cx, rownames = rownames(res$corr.X.Cx),
             colnames = colnames(res$corr.X.Cx))
    file.name = paste(output_path,'/', prefix,'_',condition, '_',postfix,'_corrYCy.mat',sep="")
    writeMat(file.name, BehCy = res$corr.Y.Cy, rownames = rownames(res$corr.Y.Cy),
             colnames = colnames(res$corr.Y.Cy))
    # write std coef and raw coef
    file.name = paste(output_path,'/', prefix,'_',condition, '_',postfix,'_coef.mat',sep="")
    writeMat(file.name, xcoef.std = res$xcoef.std, xcoef.raw = res$xcoef.raw,
             ycoef.std = res$ycoef.std, ycoef.raw = res$ycoef.raw, 
             xrownames = rownames(res$corr.X.Cx), yrownames=rownames(res$corr.Y.Cy))
    
    df_all = cbind(res$xcoef.std,res$corr.X.Cx)
    
    num_modes <- ncol(df_all) / 2
    # Generate names
    xcoef_names <- paste0("xcoef_std_m", 1:num_modes)
    xloading_names <- paste0("xloading_m", 1:num_modes)
    # Assign column names
    colnames(df_all) <- c(xcoef_names, xloading_names)
    
    file.name = paste(output_path,'/', prefix,'_',condition, '_',postfix,'_coef.csv',sep="")
    write.csv(df_all, file.name)
  }
}

