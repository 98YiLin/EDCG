### BW_Exploring Drought-response Crucial Genes in Sorghum
### Based on the P value and the expressions of genes, an example of calculating the posterior probabilities
mv=read.mv("MV_pvalue.csv") 
exp.data=read.expression("treatment.csv","control.csv")
mv.exp.data=overlay.exp.mv(exp.data,mv)
network.data=get.network(mv.exp.data)
network.bundle=dichotimize.coexpression(mv.exp.data)  
results=icmicm(network.bundle)  
write.csv(results,"Posterior probabilities.csv")
write.csv(network.bundle[[1]],"geneID.csv")
write.csv(network.bundle[[2]],"RewireMatrix.csv")
write.csv(network.bundle[[3]],"EdgeVector.csv")
write.csv(network.data[[1]],"Treatment_GDCN.csv")
write.csv(network.data[[2]],"Control_GDCN.csv")