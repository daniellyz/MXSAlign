
custom.dist <- function(my.list, my.function) {
  n <- length(my.list)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- names(my.list)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      mat[i,j] <- my.function(my.list[i],my.list[j])
    }}
  return(as.dist(mat))
}

ppm_distance<-function(x,y){
  return((x-y)/y*1000000)}

ppm_distance2<-function(x,y){
  return(abs((x-y)/y*1000000))}

cut_mass_list<-function(samples,tol){ 
  # Samples: whole data matrix
  # We align only masses
  
  samples=samples[order(samples[,3]),] # Order by mass
  N=nrow(samples)
  
  if (N>1){
  # First cut: Neighbouring masses < tol, temporary features
  f=1
  feature=c(1,rep(0,N-1)) # which mass feature it is
  
  for (k in 1:(N-1)){
    diff=(samples[k+1,3]-samples[k,3])/samples[k,3]*1000000
    if (diff<=tol){
      feature[k+1]=f}
    else{
      f=f+1
      feature[k+1]=f}}
  
  # Second cut: hierarchical clustering 
  
  feature_list=unique(feature)
  new_mass_feature=rep("0",N)
  for (f in feature_list){
    valid=which(feature==f)
    if (length(valid)>2){ 
      dis <- custom.dist(samples[valid,3], ppm_distance) # ppm distance
      discriminated<-cutree(hclust(dis,method="average"),h=tol)
      new_mass_feature[valid]=paste(f,discriminated,sep="-")}
    
    else { # If only 1-2 features no HCT is provided
      new_mass_feature[valid]= paste(f,"1",sep="-")}}}

  else{new_mass_feature="1-1"}
  
  return(cbind(samples,new_mass_feature))   
}

cut_RT_list<-function(samples,tol,t1,t2,tol2){ 
  
  # Samples: whole data matrix
  # We align RTs in each sample
  
  samples=samples[order(samples[,4]),]
  
  N=nrow(samples)
  if (N>1){
    # First cut: Neighbouring RT < tol, temporary features
    f=1
    feature=c(1,rep(0,N-1)) # which RT feature it is
    for (k in 1:(N-1)){
      diff=samples[k+1,4]-samples[k,4]
      if ((samples[k,4]>=t1) & (samples[k,4]<=t2)){tol=tol2}
      if (diff<=tol){
        feature[k+1]=f}
      else{
        f=f+1
        feature[k+1]=f}}
    
    # Second cut: hierarchical clustering 
    
    feature_list=unique(feature)
    new_RT_feature=rep("0",N)
    for (f in feature_list){
      valid=which(feature==f)
      if (length(valid)>2){ 
        dis=dist(samples[valid,4])
        discriminated<-cutree(hclust(dis,method="average"),h=tol)
        new_RT_feature[valid]=paste(f,discriminated,sep="-")}
      
      else { # If only 1-2 features no HCT is provided
        new_RT_feature[valid]= paste(f,"1",sep="-")}}}
  
  else{new_RT_feature="1-1"}
  return(cbind(samples,new_RT_feature))   # Column + which group labelled by "f1"-"f2"
}

matrix_generator<-function(samples,filenames,mtol,rtol,mode,t1,t2,tol2){
  
  # Global function
  # Find molecular features,average,evaluate mz/RT deviation and generate final matrix
  
  # Mass first:
  if (mode==1){
    masslist_cut=cut_mass_list(samples,mtol) # Cut according to mass
    masslist_cut_list=split(masslist_cut,masslist_cut[,6]) # Split into mass groups
    masslist_cut_list_cut=lapply(masslist_cut_list,function(x) cut_RT_list(x,rtol,t1,t2,tol2)) # Each mass groups splitted to RT
    samples=do.call(rbind,masslist_cut_list_cut)} # Combine back everything
  
  # RT first:
  if (mode==2){
    RTlist_cut=cut_RT_list(samples,rtol,t1,t2,tol2) # Cut according to RT
    RTlist_cut_list=split(RTlist_cut,RTlist_cut[,6]) # Split into RT groups
    RTlist_cut_list_cut=lapply(RTlist_cut_list,function(x) cut_mass_list(x,mtol)) # Each mass groups splitted to RT
    samples=do.call(rbind,RTlist_cut_list_cut)} # Combine back everything
  
  NBS=length(filenames)
  new_molecular_feature_label=paste(samples[,6],samples[,7],sep="-")
  samples_splitted=split(samples[,1:5],new_molecular_feature_label)
  N_feature=length(samples_splitted) # Number of molecular features
  I_matrix=matrix(0,N_feature,NBS) # Intensity matrix
  mass_list=rep(0,N_feature) # Averaged mass list
  RT_list=rep(0,N_feature) # Averaged RT list
  mass_dev=matrix(NA,N_feature,4) # Min + Max mass deviation, which sample
  RT_dev=matrix(NA,N_feature,4) # Min + Max RT deviation, which sample
  mass_dev_matrix=matrix(NA,N_feature,NBS) # Mass deviations in each sample
  RT_dev_matrix=matrix(NA,N_feature,NBS) # RT deviations in each sample
  
  for (i in 1:N_feature){
    sub_sample=samples_splitted[[i]]
    mass_list[i]=mean(sub_sample[,3])
    RT_list[i]=mean(sub_sample[,4])
    if (nrow(sub_sample)>1){
      ppm1=(min(sub_sample[,3])-mass_list[i])/mass_list[i]*1000000
      ppm2=(max(sub_sample[,3])-mass_list[i])/mass_list[i]*1000000
      mass_dev[i,1:2]=c(ppm1,ppm2)
      Sampleid=sub_sample[,2]
      mass_dev[i,3:4]=c(Sampleid[which.min(sub_sample[,3])],Sampleid[which.max(sub_sample[,3])])
      RT_dev[i,1:2]=c(min(sub_sample[,4])-RT_list[i],max(sub_sample[,4])-RT_list[i])
      RT_dev[i,3:4]=c(Sampleid[which.min(sub_sample[,4])],Sampleid[which.max(sub_sample[,4])])}
      for (j in 1:nrow(sub_sample)){
          I_matrix[i,sub_sample[j,2]]= I_matrix[i,sub_sample[j,2]]+sub_sample[j,5]
          mass_dev_matrix[i,sub_sample[j,2]]=sub_sample[j,3]-mass_list[i]
          RT_dev_matrix[i,sub_sample[j,2]]=sub_sample[j,4]-RT_list[i]
     }
  }
  
  colnames(I_matrix)=filenames
  matrix_generated=cbind(mass=mass_list,RT=RT_list,I_matrix)
  colnames(mass_dev)=c("Min_mz","Max_mz","Min_Sample","Max_Sample")
  colnames(RT_dev)=c("Min_RT","Max_RT","Min_Sample","Max_Sample")
  return(list(data=matrix_generated,mdev=mass_dev,rdev=RT_dev,mass_dev_matrix=mass_dev_matrix,RT_dev_matrix=RT_dev_matrix))
}



multiTitle <- function(...){ 
  ### 
  ### multi-coloured title 
  ### 
  ### examples: 
  ###  multiTitle(color="red","Traffic", 
  ###             color="orange"," light ", 
  ###             color="green","signal") 
  ### 
  ### - note triple backslashes needed for embedding quotes: 
  ### 
  ###  multiTitle(color="orange","Hello ", 
  ###             color="red"," \\\"world\\\"!") 
  ### 
  ### Barry Rowlingson <[hidden email]> 
  ### 
  l = list(...) 
  ic = names(l)=='color' 
  colors = unique(unlist(l[ic])) 
  
  for(i in colors){ 
    color=par()$col.main 
    strings=c() 
    for(il in 1:length(l)){ 
      p = l[[il]] 
      if(ic[il]){ # if this is a color: 
        if(p==i){  # if it's the current color 
          current=TRUE 
        }else{ 
          current=FALSE 
        } 
      }else{ # it's some text 
        if(current){ 
          # set as text 
          strings = c(strings,paste('"',p,'"',sep="")) 
        }else{ 
          # set as phantom 
          strings = c(strings,paste("phantom(\"",p,"\")",sep="")) 
        } 
      } 
    } # next item 
    ## now plot this color 
    prod=paste(strings,collapse="*") 
    express = paste("expression(",prod,")",sep="") 
    e=eval(parse(text=express)) 
    title(e,col.main=i) 
  } # next color 
  return() 
} 

enhanced_alignment<-function(newdata,mtol,rtol){
  
  newdata=data.frame(ID=paste0("T",1:nrow(newdata)),newdata,stringsAsFactors = F)
  ref_ID=newdata[,1]
  ref_mass=newdata[,2]
  ref_rt=newdata[,3]
  N=nrow(newdata)
  C=ncol(newdata)
  duplicated_index=c()
  duplicated=c()
  
  for (i in 1:N){
    mass=newdata[i,2]
    rt=newdata[i,3]
    dist_mass=ppm_distance2(mass,ref_mass)
    dist_rt=abs(rt-ref_rt)
    
    valid_mass=which(dist_mass<=mtol)
    valid_rt=which(dist_rt<=rtol)
    valid_both=intersect(valid_mass,valid_rt)
    
    if (length(valid_both)>1){
      duplicated_index=c(duplicated_index,valid_both)
      c1=c(ref_ID[valid_both[1]]) # Fusion 2 duplicated features
      c23=c(mean(ref_mass[valid_both]),mean(ref_rt[valid_both]))
      careas=colSums(newdata[valid_both,4:C])
      tmp=c(c1,c23,unname(careas))
      duplicated=rbind(duplicated,tmp)}}
  
  if (!is.null(duplicated)){
  rownames(duplicated)=NULL
  duplicated=duplicated[!duplicated(duplicated[,1]),]
  duplicated=data.frame(duplicated,stringsAsFactors = F)  
  duplicated[,2:C]=apply(duplicated[,2:C],2,as.numeric)
  not_duplicated=newdata[setdiff(1:N,duplicated_index),] 
  colnames(duplicated)=colnames(not_duplicated)
  rownames(duplicated)=rownames(not_duplicated)=NULL
  newdata_combined=rbind(duplicated,not_duplicated)

  valid=match(newdata_combined[,1],newdata[,1])
  newdata_combined=newdata_combined[,2:C]}
  
  else{
    valid=1:N
    newdata_combined=newdata[,2:C]
  }
  
  return(list(output=newdata_combined,valid=valid))
}


