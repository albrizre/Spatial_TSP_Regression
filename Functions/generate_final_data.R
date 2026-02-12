# Data loading

load("Results/data_clust.rda")
aux=data
load("Results/data_hom.rda")
data=rbind(aux,data)
aux=data
load("Results/data_reg.rda")
data=rbind(aux,data)
# Add name instances
data$instance_name=c()
for (i in 1:50){
  data$instance_name[i]=paste0("clust_ins_",i)
}
for (i in 51:100){
  data$instance_name[i]=paste0("hom_ins_",i-50)
}
for (i in 101:150){
  data$instance_name[i]=paste0("reg_ins_",i-100)
}
data$num=1:150

# Loading Concorde results

concorde_results=openxlsx::read.xlsx("Results/results_concorde.xlsx")
colnames(concorde_results)=c("instance_name","opt_concorde")
data=merge(data,concorde_results,by="instance_name")
data$opt_concorde=as.numeric(data$opt_concorde)
data$opt_concorde[is.na(data$opt_concorde)]=0
data$opt_concorde=data$opt_concorde/10^6
# Replace opt by opt_concorde (if non-NA, that is greater than 0)
data$opt[data$opt_concorde>0]=data$opt_concorde[data$opt_concorde>0]
# Reorder data
data=data[order(data$num),]
