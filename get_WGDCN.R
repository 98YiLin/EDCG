### BW_Exploring Drought-response Crucial Genes in Sorghum
### Constructing weighted gene differential co-expression networks(WGDCN)
#Reference1: Wang, P. and Wang, D. (2021). Gene differential co-expression networks based on RNA-seq: construction and its applications. Preprint at IEEE/ACM Trans. Comput. Biol. Bioinformat., 10.1109/TCBB.2021.3103280.
#Reference2: Hou, L., Chen, M., Zhang, C.K., Cho, J. and Zhao, H. (2014). Guilt by rewiring: gene prioritization through network rewiring in genome wide association studies. Hum. Mol. Genet. 23, 2780-2790. 10.1093/hmg/ddt668.
get.network<-function(mv.exp.data){
  treatment.exp<-mv.exp.data[[3]]   
  control.exp<-mv.exp.data[[4]]
  num_gene=nrow(treatment.exp) 
  n1<-ncol(treatment.exp)
  rho=0.5
  gr1=matrix(nrow=num_gene,ncol=num_gene)
  T1=matrix(nrow=num_gene,ncol=n1)
  for (i in 1:num_gene){
    for (j in 1:num_gene) {
      for (q in 1:n1) {
        T1[j,q]=treatment.exp[j,q]-treatment.exp[i,q]  
      }
    }
    jc11=min(min(abs(t(T1))))
    jc21=max(max(abs(t(T1))))
    ksi1=(jc11+rho*jc21)/(abs(T1)+rho*jc21)
    rt1=colMeans(t(ksi1))
    gr1[i,]=rt1
  }
  r1=(gr1+t(gr1))/2;diag(r1)=0  
  n2<-ncol(control.exp)  
  gr2=matrix(nrow=num_gene,ncol=num_gene)
  T2=matrix(nrow=num_gene,ncol=n2)
  for (i in 1:num_gene){
    for (j in 1:num_gene) {
      for (q in 1:n2) {
        T2[j,q]=control.exp[j,q]-control.exp[i,q]  
      }
    }
    jc12=min(min(abs(t(T2))))
    jc22=max(max(abs(t(T2))))
    ksi2=(jc12+rho*jc22)/(abs(T2)+rho*jc22)
    rt2=colMeans(t(ksi2))
    gr2[i,]=rt2
  }
  r2=(gr2+t(gr2))/2;diag(r2)=0  
  rewire<-abs(r1-r2)    
  g<-mv.exp.data[[2]]                
  rownames(r1)<-g; colnames(r1)<-g
  rownames(r2)<-g; colnames(r2)<-g
  rownames(rewire)<-g; colnames(rewire)<-g
  network.data<-list(r1,r2,rewire)
  network.data
}
dichotimize.coexpression<-function(mv.exp.data){
  network.data<-get.network(mv.exp.data)  
  mv<-mv.exp.data[[1]]  
  treatment.gcc<-network.data[[1]]
  control.gcc<-network.data[[2]]
  A=0.9   
  netwk<-matrix(0,ncol=ncol(treatment.gcc),nrow=nrow(treatment.gcc))  
  g<-mv.exp.data[[2]]                                
  rownames(netwk)<-g                                   
  colnames(netwk)<-g                                   
  seclet=matrix(ncol=ncol(treatment.gcc),nrow=nrow(treatment.gcc))
  seclet=(abs(treatment.gcc)>=A)+(abs(control.gcc)>=A)
  netwk[which(seclet==1)]<-1
  degree<-rowSums(netwk)                              
  gmrf<-g[degree>0]                                   
  rewire<-network.data[[3]][gmrf,gmrf]                 
  netwk<-netwk[gmrf,gmrf]                             
  netwk1<-netwk
  netwk1[lower.tri(netwk1)]<-0                         
  nNodes<-length(gmrf)
  genesymbol<-c(1:nNodes)
  mv<-mv[gmrf,]
  edges <- expand.grid(genesymbol,genesymbol)[as.vector(netwk1)==1,]
  NEI<-function(ii){which(netwk[ii,]==1)}  
  neighbors<-lapply(1:nNodes,"NEI")        
  names(neighbors)<-as.character(1:nNodes)
  edge.weight<-as.vector(rewire)[as.vector(netwk1)==1]
  edges<-cbind(edges,edge.weight)          
  network.bundle<-list(gmrf,rewire,edges,neighbors,mv)
  names(network.bundle)<-c("GeneID","RewireMatrix","EdgeVector","NeighborList","MV")  
  network.bundle                                     
}  
