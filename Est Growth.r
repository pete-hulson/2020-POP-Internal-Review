
#====================================================================================================
#================ Get Working Directories, define format vectors, define data call parameters
#====================================================================================================

######### Get directories

path<-getwd()
path_DAT<-paste0(path,"/Data")
path_MDL<-paste0(path,"/Model")
path_RES<-paste0(path,"/Results")
path_RES_Ann<-paste0(path,"/Results/Annual Results")

source(paste0(path,"/Growth Fcns.R"))

#====================================================================================================
#===== Read in data
#====================================================================================================

lfreq<-read.csv(paste0(path_DAT,"/lfreq_full.csv"))
specimen<-read.csv(paste0(path_DAT,"/specimen_full.csv"))
strata<-read.csv(paste0(path_DAT,"/strata.csv"))

# Subset to data after 1990
lfreq<-subset(lfreq,lfreq$YEAR>=1990)
specimen<-subset(specimen,specimen$YEAR>=1990)

# Remove NA's from age data
specimen<-subset(specimen,is.na(specimen$AGE)==FALSE)

# Get regional data
str_WGOA<-strata$STRATUM[which(strata$REGULATORY_AREA_NAME == "WESTERN GOA")]
str_CGOA<-strata$STRATUM[which(strata$REGULATORY_AREA_NAME == "CENTRAL GOA")]
str_EGOA<-strata$STRATUM[which(strata$REGULATORY_AREA_NAME == "EASTERN GOA")]

SUBREGION<-vector("character",length = length(lfreq$STRATUM))
SUBREGION[which(lfreq$STRATUM %in% str_WGOA)]<-"WGOA"
SUBREGION[which(lfreq$STRATUM %in% str_CGOA)]<-"CGOA"
SUBREGION[which(lfreq$STRATUM %in% str_EGOA)]<-"EGOA"
lfreq<-cbind(lfreq,SUBREGION)
lfreq_WGOA<-subset(lfreq,lfreq$SUBREGION == "WGOA")
lfreq_CGOA<-subset(lfreq,lfreq$SUBREGION == "CGOA")
lfreq_EGOA<-subset(lfreq,lfreq$SUBREGION == "EGOA")

SUBREGION<-vector("character",length = length(specimen$STRATUM))
SUBREGION[which(specimen$STRATUM %in% str_WGOA)]<-"WGOA"
SUBREGION[which(specimen$STRATUM %in% str_CGOA)]<-"CGOA"
SUBREGION[which(specimen$STRATUM %in% str_EGOA)]<-"EGOA"
specimen<-cbind(specimen,SUBREGION)
specimen_WGOA<-subset(specimen,specimen$SUBREGION == "WGOA")
specimen_CGOA<-subset(specimen,specimen$SUBREGION == "CGOA")
specimen_EGOA<-subset(specimen,specimen$SUBREGION == "EGOA")


#====================================================================================================
#===== Model Input Parameters/Results tables
#====================================================================================================

yrs<-sort(unique(specimen$YEAR))
ages_M<-seq(2,35)
nages_M<-length(ages_M)

W_ests_GOA<-matrix(nrow = nages_M+4,ncol = length(yrs)+1)
rownames(W_ests_GOA)<-c("Winf","k","t0","beta",ages_M)
colnames(W_ests_GOA)<-c("All",yrs)

W_ests_WGOA<-matrix(nrow = nages_M+4,ncol = length(yrs)+1)
rownames(W_ests_WGOA)<-c("Winf","k","t0","beta",ages_M)
colnames(W_ests_WGOA)<-c("All",yrs)

W_ests_CGOA<-matrix(nrow = nages_M+4,ncol = length(yrs)+1)
rownames(W_ests_CGOA)<-c("Winf","k","t0","beta",ages_M)
colnames(W_ests_CGOA)<-c("All",yrs)

W_ests_EGOA<-matrix(nrow = nages_M+4,ncol = length(yrs)+1)
rownames(W_ests_EGOA)<-c("Winf","k","t0","beta",ages_M)
colnames(W_ests_EGOA)<-c("All",yrs)


L_ests_GOA<-matrix(nrow = nages_M+3,ncol = length(yrs)+1)
rownames(L_ests_GOA)<-c("Linf","k","t0",ages_M)
colnames(L_ests_GOA)<-c("All",yrs)

L_ests_WGOA<-matrix(nrow = nages_M+3,ncol = length(yrs)+1)
rownames(L_ests_WGOA)<-c("Linf","k","t0",ages_M)
colnames(L_ests_WGOA)<-c("All",yrs)

L_ests_CGOA<-matrix(nrow = nages_M+3,ncol = length(yrs)+1)
rownames(L_ests_CGOA)<-c("Linf","k","t0",ages_M)
colnames(L_ests_CGOA)<-c("All",yrs)

L_ests_EGOA<-matrix(nrow = nages_M+3,ncol = length(yrs)+1)
rownames(L_ests_EGOA)<-c("Linf","k","t0",ages_M)
colnames(L_ests_EGOA)<-c("All",yrs)



#====================================================================================================
#================ Run functions to estimate growth - GOA wide
#====================================================================================================

# Combined years

wt<-get_waa(specimen,lfreq,"AllY",nages_M,ages_M,"AllR",path_RES_Ann)
W_ests_GOA[,1]<-wt

ln<-get_szaa(specimen,lfreq,"AllY",nages_M,ages_M,"AllR",path_RES_Ann)
L_ests_GOA[,1]<-ln

# Annual

for (y in 1:length(yrs)){
  
  spec_y<-subset(specimen,specimen$YEAR == yrs[y])
  lfr_y<-subset(lfreq,lfreq$YEAR == yrs[y])
  
  wt<-get_waa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"AllR",path_RES_Ann)
  W_ests_GOA[,y+1]<-wt
  
  ln<-get_szaa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"AllR",path_RES_Ann)
  L_ests_GOA[,y+1]<-ln
  
  print(y)
  
}

write.csv(W_ests_GOA,paste0(path_RES,"/Wt_ests_GOA.csv"))
write.csv(L_ests_GOA,paste0(path_RES,"/Lnth_ests_GOA.csv"))

#====================================================================================================
#================ Run functions to estimate growth - by region
#====================================================================================================

### WGOA

# Combined years

wt<-get_waa(specimen_WGOA,lfreq_WGOA,"AllY",nages_M,ages_M,"WGOA",path_RES_Ann)
W_ests_WGOA[,1]<-wt

ln<-get_szaa(specimen_WGOA,lfreq_WGOA,"AllY",nages_M,ages_M,"WGOA",path_RES_Ann)
L_ests_WGOA[,1]<-ln

# Annual

for (y in 1:length(yrs)){
  
  spec_y<-subset(specimen_WGOA,specimen_WGOA$YEAR == yrs[y])
  lfr_y<-subset(lfreq_WGOA,lfreq_WGOA$YEAR == yrs[y])
  
  wt<-get_waa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"WGOA",path_RES_Ann)
  W_ests_WGOA[,y+1]<-wt
  
  ln<-get_szaa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"WGOA",path_RES_Ann)
  L_ests_WGOA[,y+1]<-ln
  
}

write.csv(W_ests_WGOA,paste0(path_RES,"/Wt_ests_WGOA.csv"))
write.csv(L_ests_WGOA,paste0(path_RES,"/Lnth_ests_WGOA.csv"))




### CGOA

# Combined years

wt<-get_waa(specimen_CGOA,lfreq_CGOA,"AllY",nages_M,ages_M,"CGOA",path_RES_Ann)
W_ests_CGOA[,1]<-wt

ln<-get_szaa(specimen_CGOA,lfreq_CGOA,"AllY",nages_M,ages_M,"CGOA",path_RES_Ann)
L_ests_CGOA[,1]<-ln

# Annual

for (y in 1:length(yrs)){
  
  spec_y<-subset(specimen_CGOA,specimen_CGOA$YEAR == yrs[y])
  lfr_y<-subset(lfreq_CGOA,lfreq_CGOA$YEAR == yrs[y])
  
  wt<-get_waa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"CGOA",path_RES_Ann)
  W_ests_CGOA[,y+1]<-wt
  
  ln<-get_szaa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"CGOA",path_RES_Ann)
  L_ests_CGOA[,y+1]<-ln
  
}

write.csv(W_ests_CGOA,paste0(path_RES,"/Wt_ests_CGOA.csv"))
write.csv(L_ests_CGOA,paste0(path_RES,"/Lnth_ests_CGOA.csv"))


### EGOA

# Combined years

wt<-get_waa(specimen_EGOA,lfreq_EGOA,"AllY",nages_M,ages_M,"EGOA",path_RES_Ann)
W_ests_EGOA[,1]<-wt

ln<-get_szaa(specimen_EGOA,lfreq_EGOA,"AllY",nages_M,ages_M,"EGOA",path_RES_Ann)
L_ests_EGOA[,1]<-ln

# Annual

for (y in 1:length(yrs)){
  
  spec_y<-subset(specimen_EGOA,specimen_EGOA$YEAR == yrs[y])
  lfr_y<-subset(lfreq_EGOA,lfreq_EGOA$YEAR == yrs[y])
  
  if(length(spec_y$WEIGHT)>length(which(is.na(spec_y$WEIGHT)))){
    wt<-get_waa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"EGOA",path_RES_Ann)
    W_ests_EGOA[,y+1]<-wt}
  
  ln<-get_szaa(spec_y,lfr_y,yrs[y],nages_M,ages_M,"EGOA",path_RES_Ann)
  L_ests_EGOA[,y+1]<-ln
  
}

write.csv(W_ests_EGOA,paste0(path_RES,"/Wt_ests_EGOA.csv"))
write.csv(L_ests_EGOA,paste0(path_RES,"/Lnth_ests_EGOA.csv"))












