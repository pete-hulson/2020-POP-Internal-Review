

#~~~~~~~~~~~~~~~~~~~ Function to compute current mean weight
get_waa<-function(age_data,length_data,year,nages_M,ages_M,region,path_RES){
  
  
  # Get parameters
  ages<-sort(unique(age_data$AGE))
  nages<-length(ages)
  lengths<-sort(unique(age_data$LENGTH))
  nlengths<-length(lengths)
  
  # Subset out age data with no weights
  r<-which(is.na(age_data_1$WEIGHT)==TRUE)
  if(length(r)>0)
    age_data_1 <- age_data_1[-r,]
  ages<-sort(unique(age_data_1$AGE))
  nages<-length(ages)
  lengths<-sort(unique(age_data_1$LENGTH))
  nlengths<-length(lengths)
    
  # Subset to ages with >1 obs
  n_a<-table(age_data$AGE)
  r<-which(n_a<2)
  if(length(r)>0)
    n_a<-n_a[-r]
  ages<-as.numeric(names(n_a))
  nages<-length(ages)
  age_data_1 <- subset(age_data,age_data$AGE %in% ages)
  
  # Get Age-length key together
  n_al<-table(age_data_1$AGE,age_data_1$LENGTH)
  n_l<-colSums(n_al)
  r<-which(n_l<2)
  if(length(r)>0){
    n_l<-n_l[-r]
    n_al<-n_al[,-r]}
  lengths<-as.numeric(names(n_l))
  nlengths<-length(lengths)
  len_data_sub<-subset(length_data,length_data$LENGTH %in% lengths)
  N_l<-tapply(len_data_sub$FREQUENCY,len_data_sub$LENGTH,sum)
  N_al<-sweep(sweep(n_al,2,n_l,'/'),2,N_l,'*')
  
  # Get mean/variance in weight age-length key
  Wbar_la<-tapply(age_data_1$WEIGHT,list(age_data_1$AGE,age_data_1$LENGTH),FUN=mean,na.rm=TRUE)
  Wbar_la<-Wbar_la[,which(as.numeric(colnames(Wbar_la)) %in% lengths)]
  V_Wbar_la<-tapply(age_data_1$WEIGHT,list(age_data_1$AGE,age_data_1$LENGTH),FUN=var,na.rm=TRUE)/tapply(age_data_1$WEIGHT,list(age_data_1$AGE,age_data_1$LENGTH),FUN=length)
  V_Wbar_la<-V_Wbar_la[,which(as.numeric(colnames(V_Wbar_la)) %in% lengths)]
  alpha_l<-N_l/sum(N_l)
  theta_la<-sweep(n_al,2,colSums(n_al),'/')
  r_la<-sweep(theta_la,2,alpha_l,'*')
  theta_a<-rowSums(r_la)
  L<-sum(N_l)
  A_l<-colSums(n_al)
  V_r_la<-matrix(NA,nrow=nages,ncol=nlengths)
  rownames(V_r_la)<-ages
  colnames(V_r_la)<-lengths
  for(a in 1:nages){
    for(l in 1:nlengths){
      V_r_la[a,l]<-alpha_l[l]^2*theta_la[a,l]*(1-theta_la[a,l])/(A_l[l]-1)+alpha_l[l]*(theta_la[a,l]-theta_a[a])^2/L}}
  
  # Get/compile weight-at-age statistics
  Age<-ages
  SS<-tapply(age_data_1$WEIGHT,age_data_1$AGE,FUN=length)
  Wbar<-rowSums(r_la*Wbar_la,na.rm=TRUE)/rowSums(r_la)
  SD_Wbar<-vector(length=nages)  
  for(a in 1:nages){
    SD_Wbar[a]<-sqrt(sum(r_la[a,]^2*V_Wbar_la[a,]+(Wbar_la[a,]-Wbar[a])^2*V_r_la[a,],na.rm=TRUE)/theta_a[a]^2)*sqrt(length(subset(age_data_1$WEIGHT,age_data_1$AGE==ages[a])))}
  WaA_stats<-as.data.frame(cbind(Age,SS,Wbar,SD_Wbar))
  r<-which(WaA_stats$SD_Wbar==0)
  if(length(r)>0)
    WaA_stats<-WaA_stats[-r,]
  r<-which(is.na(WaA_stats$SD_Wbar)==TRUE)
  if(length(r)>0)
    WaA_stats<-WaA_stats[-r,]
  
  # Write mean weight-at-age
  write.csv(WaA_stats,paste0(path_RES,"/mean_WaA_",year,"_",region,".csv"))
  
  # Get/compile weight-at-length statistics
  lw_data<-cbind(age_data$LENGTH,age_data$WEIGHT)
  r<-which(is.na(lw_data[,2])==TRUE)
  if(length(r)>0)
    lw_data<-lw_data[-r,]
  lengths<-sort(unique(lw_data[,1]))
  lw_mdl_data<-matrix(nrow=length(lengths),ncol=3)
  colnames(lw_mdl_data)<-c("Length","Wbar","SD_Wbar")
  lw_mdl_data[,1]<-lengths
  lw_mdl_data[,2]<-tapply(lw_data[,2],lw_data[,1],mean)
  lw_mdl_data[,3]<-tapply(lw_data[,2],lw_data[,1],sd)
  r<-which(is.na(lw_mdl_data[,3])==TRUE)
  if(length(r)>0)
    lw_mdl_data<-lw_mdl_data[-r,]
  r<-which(lw_mdl_data[,3]==0)
  if(length(r)>0)
    lw_mdl_data<-lw_mdl_data[-r,]
  
  # Write mean weight-at-length
  write.csv(lw_mdl_data,paste0(path_RES,"/mean_WaL_",year,"_",region,".csv"))
  
  # Run allometric model
  DAT<-c("# Data file for allometric model of mean weight by length",
         "# Number of lengths (nlengths)",
         length(lw_mdl_data[,1]),
         "# Lengths with observed mean weight (lengths)",
         paste(lw_mdl_data[,1],collapse=" "),
         "# Observed mean weight (Wbar_obs)",
         paste(lw_mdl_data[,2],collapse=" "),
         "# SD in Observed mean weight (SD_Wbar)",
         paste(lw_mdl_data[,3],collapse=" "))
  write.table(DAT,file=paste0(path_MDL,"/Wbar@A/ALLO.DAT"),quote=F,row.names=F,col.names=F)
  setwd(paste0(path_MDL,"/Wbar@A"))
  shell("ALLO.EXE")
  PAR<-readLines(paste0(path_MDL,"/Wbar@A/allo.par"),warn=FALSE)
  beta_lw<-as.numeric(strsplit(PAR[grep("beta",PAR)+1]," ")[[1]])
  
  # Run LVBmodel and estimate mean weight
  PIN<-c("# Parameter starting values for LVB model of mean weight","# Winf","800","# k","0.1","# t0","0","# beta",as.character(beta_lw))
  write.table(PIN,file=paste0(path_MDL,"/Wbar@A/LVB.PIN"),quote=F,row.names=F,col.names=F)
  DAT<-c("# Data file for LVB model of mean weight",
         "# Number of ages (nages)",
         length(WaA_stats$Age),
         "# Ages with observed mean weight (ages)",
         paste(WaA_stats$Age,collapse=" "),
         "# Observed mean weight (Wbar_obs)",
         paste(WaA_stats$Wbar,collapse=" "),
         "# SD in Observed mean weight (Wbar_obs)",
         paste(WaA_stats$SD_Wbar,collapse=" "))
  write.table(DAT,file=paste0(path_MDL,"/Wbar@A/LVB.DAT"),quote=F,row.names=F,col.names=F)
  shell("LVB.EXE")
  REP<-readLines(paste(path_MDL,"/Wbar@A/lvb.rep",sep=""),warn=FALSE)
  Winf<-as.numeric(strsplit(REP[grep("Winf",REP)[1]]," ")[[1]][2])
  k<-as.numeric(strsplit(REP[grep("k",REP)[1]]," ")[[1]][2])
  t0<-as.numeric(strsplit(REP[grep("t0",REP)[1]]," ")[[1]][2])
  Wbar<-Winf*(1-exp(-k*(ages_M-t0)))^beta_lw
  Wbar[nages_M]<-0.5*(Wbar[nages_M]+Winf)
  Wbar<-round(Wbar,digits=1)
  Wbar_params<-cbind(Winf,k,t0,beta_lw)
  write.csv(Wbar_params,paste0(path_RES,"/Wbar_params_",year,"_",region,".csv"))
  
  c(Wbar_params,Wbar)}


#~~~~~~~~~~~~~~~~~~~ Function to compute current size-age transition matrix
get_szaa<-function(age_data,length_data,year,nages_M,ages_M,region,path_RES){
  
  ### Compute observed length-at-age statistics
  # Get parameters
  ages<-sort(unique(age_data$AGE))
  nages<-length(ages)
  lengths<-sort(unique(age_data$LENGTH))
  nlengths<-length(lengths)
  
  # Subset to ages with >1 obs
  n_a<-table(age_data$AGE)
  r<-which(n_a<2)
  if(length(r)>0)
    n_a<-n_a[-r]
  ages<-as.numeric(names(n_a))
  nages<-length(ages)
  age_data_1 <- subset(age_data,age_data$AGE %in% ages)
  lengths<-sort(unique(age_data_1$LENGTH))
  nlengths<-length(lengths)
  
  # Subset out age data with no lengths
  r<-which(is.na(age_data_1$LENGTH)==TRUE)
  if(length(r)>0)
    age_data_1 <- age_data_1[-r,]
  ages<-sort(unique(age_data_1$AGE))
  nages<-length(ages)
  lengths<-sort(unique(age_data_1$LENGTH))
  nlengths<-length(lengths)
  
  # Get Age-length key together
  n_al<-table(age_data_1$AGE,age_data_1$LENGTH)
  n_l<-colSums(n_al)
  r<-which(n_l<2)
  if(length(r)>0){
    n_l<-n_l[-r]
    n_al<-n_al[,-r]}
  lengths<-as.numeric(names(n_l))
  nlengths<-length(lengths)
  len_data_sub<-subset(length_data,length_data$LENGTH %in% lengths)
  N_l<-tapply(len_data_sub$FREQUENCY,len_data_sub$LENGTH,sum)
  N_al<-sweep(sweep(n_al,2,n_l,'/'),2,N_l,'*')

  # Get/compile observed length-at-age statistics
  Age<-ages
  SS<-tapply(age_data_1$LENGTH,age_data_1$AGE,FUN=length)
  Lbar<-rowSums(sweep(N_al,2,lengths/10,'*'))/rowSums(N_al)
  SD_Lbar<-vector(length=nages)
  for(a in 1:nages){
    SD_Lbar[a]<-sqrt(1/(sum(N_al[a,])-1)*sum(N_al[a,]*(lengths/10-Lbar[a])^2))}
  LaA_stats<-as.data.frame(cbind(Age,SS,Lbar,SD_Lbar))
  r<-which(LaA_stats$SD_Lbar<=0.01)
  if(length(r)>0)
    LaA_stats<-LaA_stats[-r,]
  r<-which(is.na(LaA_stats$SD_Lbar)==TRUE)
  if(length(r)>0)
    LaA_stats<-LaA_stats[-r,]
  
  # Write data
  write.csv(LaA_stats,paste0(path_RES,"/mean_LaA_",year,"_",region,".csv",sep=""))
  
  ### Estimate length-at-age statistics
  # Estimate mean length
  DAT<-c("# Data file for LVB model of mean length",
         "# Number of ages (nages)",
         length(LaA_stats$Age),
         "# Ages with observed mean length (ages)",
         paste(LaA_stats$Age,collapse=" "),
         "# Observed mean length (Lbar_obs)",
         paste(LaA_stats$Lbar,collapse=" "),
         "# SD in Observed mean length (Lbar_obs)",
         paste(LaA_stats$SD_Lbar,collapse=" "))
  write.table(DAT,file=paste0(path_MDL,"/Lbar@A/LVB.DAT"),quote=F,row.names=F,col.names=F)
  setwd(paste(path_MDL,"/Lbar@A",sep=""))
  shell("LVB.EXE")
  REP_Lbar<-readLines(paste(path_MDL,"/Lbar@A/LVB.REP",sep=""),warn=FALSE)
  Linf<-as.numeric(strsplit(REP_Lbar[grep("Linf",REP_Lbar)[1]]," ")[[1]][2])
  k<-as.numeric(strsplit(REP_Lbar[grep("k",REP_Lbar)[1]]," ")[[1]][2])
  t0<-as.numeric(strsplit(REP_Lbar[grep("t0",REP_Lbar)[1]]," ")[[1]][2])
  Lbar_params<-cbind(Linf,k,t0)
  write.csv(Lbar_params,paste0(path_RES,"/Lbar_params_",year,"_",region,".csv",sep=""))
  Lbar<-Linf*(1-exp(-k*(ages_M-t0)))
  Lbar[nages_M]<-0.5*(Lbar[nages_M]+Linf)

  c(Lbar_params,Lbar)}





