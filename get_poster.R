### BW_Exploring Drought-response Crucial Genes in Sorghum
### Calculating the posterior probabilities for candidate genes
#Reference: Hou, L., Chen, M., Zhang, C.K., Cho, J. and Zhao, H. (2014). Guilt by rewiring: gene prioritization through network rewiring in genome wide association studies. Hum. Mol. Genet. 23, 2780-2790. 10.1093/hmg/ddt668.
estimate.h1<-function(netwk, priorPara1=c(0.01,0.01),q1=0.9,delta=0.95){
  geneid<-netwk[[1]]
  mv<-netwk[[5]]
  tau11<-priorPara1[1]   
  tau22<-priorPara1[2]   
  neighbors<-netwk[[4]]   
  nNodes<-length(geneid)
  genesymbol<-1:nNodes
  nodesLabel<-ifelse(mv[,1]<=0.005,1,-1)    
  pdiff<-netwk[[2]]        
  dpot<-lapply(1:nNodes,function(nodei)
  { i.neighbor <- neighbors[[nodei]]
  i.neighbor.1<-is.element(genesymbol,i.neighbor) & nodesLabel==1    
  i.neighbor.2<-is.element(genesymbol,i.neighbor) & nodesLabel==-1 & pdiff[nodei,]>delta
  logitp<-tau11*sum(pdiff[nodei,]*i.neighbor.1)+tau22*sum(pdiff[nodei,]*i.neighbor.2)
  logitp   
  }
  )
  dpot<-unlist(dpot)  
  h1<-as.numeric(quantile(dpot,q1)*(-1))  
  h1  
}
MRF.condProb <- function(nodei,nodesLabel,netwk,priorPara1,h1,delta=0.95)
{
  neighbors<-netwk[[4]]
  genesymbol<-1:length(netwk[[1]])
  pdiff<-netwk[[2]]
  tau11<-priorPara1[1]  
  tau22<-priorPara1[2]  
  i.neighbor <- neighbors[[nodei]]                   
  i.neighbor.1<-is.element(genesymbol,i.neighbor) & nodesLabel==1
  i.neighbor.2<-is.element(genesymbol,i.neighbor) & nodesLabel==-1 & pdiff[nodei,]>delta
  logitp<-(h1+tau11*sum(pdiff[nodei,]*i.neighbor.1)+tau22*sum(pdiff[nodei,]*i.neighbor.2))  
  1/(1+exp(-logitp))     
}
updateP<-function(netwk,nodesLabel,priorPara1=c(0.01,0.01),h1)
{
  tau11<-priorPara1[1]   
  tau22<-priorPara1[2]   
  nNodes<-length(netwk[[1]])
  likRatio<-netwk[[5]][,2]
  prob<-numeric(nNodes)  
  for(nodeInd in 1:nNodes )
  {
    condProb<-MRF.condProb(nodeInd,nodesLabel,netwk,priorPara1,h1) 
    prob[nodeInd]<-(likRatio[nodeInd]*condProb)/(likRatio[nodeInd]*condProb+(1-condProb))    
  }
  prob
}
graphPrior <- function(netwk,nodesLabel,priorPara1=c(0.01,0.01),delta=0.95,h1)
{
  tau11<-priorPara1[1] 
  tau22<-priorPara1[2] 
  edges<-netwk[[3]]
  edge11 <- nodesLabel[edges[,1]]==1 & nodesLabel[edges[,2]]==1  
  edge22<-nodesLabel[edges[,1]]==-1 & nodesLabel[edges[,2]]==-1 & edges[,3]>delta  
  U<-h1*sum(nodesLabel==1)-tau22*sum(edges[edge22,3])    
  U
}
ichip.llk <- function(netwk,nodesLabel,priorPara2=c(3,1,10,1))
{
  y<-netwk[[5]][,2]           
  mubar <- priorPara2[1]            
  a <- priorPara2[2]                
  nu <- priorPara2[3]              
  lambda <- priorPara2[4]          
  scale0 <- sqrt(lambda*(a+1)/a)    
  lk<-ifelse(nodesLabel==1,dt((y-mubar)/scale0,nu)/scale0,dnorm(y,0,1) )   
  sum(log(lk))   
}
pseudo.llk<-function(netwk,nodesLabel,priorPara1=c(0.01,0.01),priorPara2=c(3,1,10,1),h1)
{
  ichip.llk(netwk,nodesLabel)+graphPrior(netwk,nodesLabel,priorPara1,0.95,h1)
}       
icmicm<-function(netwk,priorPara1=c(0.01,0.01),priorPara2=c(3,1,10,1),delta=0.95,q1=0.9)
{
  h1<-estimate.h1(netwk)  
  nNodes<-length(netwk[[1]])
  initLabel<-NULL         
  currentNodesLabel <- if(is.null(initLabel)) {2*rbinom(nNodes,1,0.5)-1}  
  likRatio<-netwk[[5]][,2]                                               
  prevNodesLabel<-currentNodesLabel
  iteNum=1;iteT=100;converged=F   
  randsearch <- function(initNode,currentNodesLabel)  
  {
    temp.label<-currentNodesLabel   
    for(nodei in sample(nNodes))
    {
      condProb<-MRF.condProb(nodei,temp.label,netwk,priorPara1,h1)  
      prob<-(likRatio[nodei]*condProb)/(likRatio[nodei]*condProb+(1-condProb))
      temp.label[nodei]<-sign(prob-0.5)  
    }
    temp.label
  }
  while (iteNum<iteT & !converged)   
  {
    prevNodesLabel<-currentNodesLabel
    currentNodesLabel<-randsearch(sample(1:nNodes,1),currentNodesLabel)
    if(all(currentNodesLabel==prevNodesLabel)) 
    { converged<-T }
    a<-pseudo.llk(netwk,prevNodesLabel,priorPara1,priorPara2,h1)   
    b<-pseudo.llk(netwk,currentNodesLabel,priorPara1,priorPara2,h1)
    if(a>b) 
    { currentNodesLabel<-prevNodesLabel
    converged<-T
    }
    iteNum <- iteNum+1
  }
  post.mv<-updateP(netwk,currentNodesLabel,priorPara1,h1)
  outout<-data.frame(netwk[[1]],netwk[[5]][,1],post.mv)
  colnames(outout)<-c("GeneID","MV","Posterior")
  outout  
}
