### BW_Exploring Drought-response Crucial Genes in Sorghum
### Reading the P value and the expression of genes
#Reference: Hou, L., Chen, M., Zhang, C.K., Cho, J. and Zhao, H. (2014). Guilt by rewiring: gene prioritization through network rewiring in genome wide association studies. Hum. Mol. Genet. 23, 2780-2790. 10.1093/hmg/ddt668.
pval.class <- function(pv) {
  pv[pv<=.Machine$double.eps] <- .Machine$double.eps  
  nodeD <- pv                                        
  class(nodeD) <- "pval"          
  attr(nodeD,'p') <- pv             
  attr(nodeD,'y') <- qnorm(1-pv);  
  nodeD
}  
likpw.pval <- function( nodesData, nodesLabel, priorPara2=NULL) {
  y <- attr(nodesData,'y') 
  mubar <- priorPara2[1]
  a <- priorPara2[2]
  nu <- priorPara2[3]
  lambda <- priorPara2[4]
  scale0 <- sqrt(lambda*(a+1)/a)
  ifelse(nodesLabel==1, dt((y-mubar)/scale0,nu)/scale0, dnorm(y,0,1))
}   
read.mv<-function(mv.name,priorPara2=c(3,1,10,1)){
  pmv<-read.csv(file=mv.name,row.names=1,header=T)
  gname<-rownames(pmv)
  p<-pmv[,1]    
  names(p)<-gname 
  p=p[p<=0.01] 
  gname<-names(p)
  obsAsso<-pval.class(p)  
  NN<-length(p)
  lik.1<-likpw.pval(obsAsso, rep(1,NN), priorPara2)   
  lik.0<-likpw.pval(obsAsso, rep(-1,NN))              
  likRatio<-lik.1/lik.0                                
  mv<-data.frame(p,likRatio)
  names(mv)<-c("pvalue","LLR")
  rownames(mv)<-gname
  return(mv)
}
read.expression<-function(treatment.name,control.name)
{
  treatment.exp<-read.csv(file=treatment.name,header=T,row.names=1)
  control.exp<-read.csv(file=control.name,header=T,row.names=1)
  g.microarray<-rownames(treatment.exp)
  exp.meta<-list(g.microarray,treatment.exp,control.exp)
}
overlay.exp.mv<-function(exp.meta,mv)
{
  g.microarray<-exp.meta[[1]]  
  g.mv<-rownames(mv)
  g<-intersect(g.microarray,g.mv)  
  mv<-mv[g,]                     
  treatment.exp<-exp.meta[[2]][g,]     
  control.exp<-exp.meta[[3]][g,]    
  mv.exp.meta<-list(mv,g,treatment.exp,control.exp)
  return(mv.exp.meta)
}