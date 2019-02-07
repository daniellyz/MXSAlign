library(shiny)
library(V8)
library(d3heatmap)
library(shinyjs)
require(DT, quietly = TRUE) 
source("scripts.r")

# Define server logic required to draw a histogram
shinyServer(function(input, output,clientData, session) {

## Load data transform and check format:
  
load_files <-eventReactive(input$goButton,{
  
  myfiles=list()
  samples=NULL
  mms=""
  NBS=length(input$file1[,1]) # Number of samples
  filenames=input$file1[['name']]
  
# Check each file and report error 
  if ((NBS<2) || (is.null(NBS))){mms="At least 2 files needed for alignment!"}
  else {
  for(i in 1:NBS){
    name=filenames[i]
    file1 <- read.table(input$file1$datapath[i],sep=",",dec=".",header=F,stringsAsFactors=F)
    if (ncol(file1)!=5){mms=paste0("The file ",name," must contain 5 columns!")}
    else if (!is.numeric(as.matrix(file1))){mms=paste0("Numeric matrix is needed in file ",name)}
    else{
      if (input$Mass=="Max. m/z"){
        myfiles[[i]]=file1[,c(1,5,3,4)]}
      if (input$Mass=="Molecular weight"){
        if (input$Scan=="Positive ion mode"){file1[,2]=file1[,2]+1.007276} # Correction
        if (input$Scan=="Negative ion mode"){file1[,2]=file1[,2]-1.007276}
        myfiles[[i]]=file1[,c(1,2,3,4)]
      }
      ### Remove all rows that contain missing values
      myfiles[[i]]=na.omit(myfiles[[i]])
    }
  }}
  
  # Combine file to sample table:  
 if (mms==""){
    for (i in 1:NBS){
      newdata=myfiles[[i]]
      newdata[,1]=i
      samples=rbind(samples,newdata)
    }
  colnames(samples)=c("Sample_id","Mass","RT","Area")
  samples=samples[samples[,3]>input$Calibrant & samples[,3]<input$End,]
  samples=cbind(ID=1:nrow(samples),samples)
  mms="Data format valid!"}
  
  list(samples=samples,mms=mms,filenames=filenames,NBS=NBS)
})


## Execute alignment or tune:

alignment <-eventReactive(input$goButton,{

  final=NULL
  tune_nbf=NULL
  tune_list=list()
  loss=NULL
  mtol_test=rtol_test=NULL
  
  samples=load_files()$samples
  filenames=load_files()$filenames
  
  if (!is.null(samples)){
     if (input$Fixed){
       t1=input$range_tt[1]
       t2=input$range_tt[2]
       tol2=input$tol2}
     else {
       t1=input$Calibrant
       t2=input$End
       tol2=input$rtol}
     if (input$Mode=="Fixed tolerance window"){
       if (input$Order=="Mass first"){
           final=try(matrix_generator(samples,filenames,input$mtol,input$rtol,1,t1,t2,tol2),silent=T)}
       if (input$Order=="Retention time first"){
           final=try(matrix_generator(samples,filenames,input$mtol,input$rtol,2,t1,t2,tol2),silent=T)}
       
       if (input$enhance){
          n0=nrow(final$data)
          tmp=enhanced_alignment(final$data,input$mtol,tol2)
          final$data=tmp$output
          final$mdev=final$mdev[tmp$valid,]
          final$rdev=final$rdev[tmp$valid,]
          final$mass_dev_matrix=final$mass_dev_matrix[tmp$valid,]
          final$RT_dev_matrix=final$RT_dev_matrix[tmp$valid,]
          loss=n0-nrow(final$data) # Loss by enhanced alignment
          }
       }
     
     if (input$Mode=="Parameter Tuning (Time-consuming)"){
       mtol_test=seq(input$range_mtol[1],input$range_mtol[2],input$seg_mtol)
       rtol_test=seq(input$range_rtol[1],input$range_rtol[2],input$seg_rtol)
       M=length(mtol_test)
       R=length(rtol_test)
       tune_nbf=matrix(0,M,R)
       for (i in 1:M){
         for (j in 1:R){
           if (input$Order=="Mass first"){
             if (!input$Fixed){tol2=rtol_test[j]}
             final2=try(matrix_generator(samples,filenames,mtol_test[i],rtol_test[j],1,t1,t2,tol2),silent=T)}
           if (input$Order=="Retention time first"){
             final2=try(matrix_generator(samples,filenames,mtol_test[i],rtol_test[j],2,t1,t2,tol2),silent=T)}
           tune_nbf[i,j]=nrow(final2$data) # Save number of features detected
           tune_list[[paste0(mtol_test[i],"-",rtol_test[j])]]=final2}
         }}
       }
  if (class(final)!="try-error" & !is.null(final)){mms="Alignment succeeded! Please go to tab-panel C) to check the results!"}
  else if (!is.null(tune_nbf)){mms="Tuning succeeded! Please go to tab-panel B) to check the results!"}
  else {mms="Alignment failed!"}
  list(final=final,mms=mms,tune_nbf=tune_nbf,tune_list=tune_list,mtol_test=mtol_test,rtol_test=rtol_test,loss=loss)
})

  
# Execute & Deliver message 

observeEvent(input$goButton,{
  
  withProgress({
    setProgress(message="Check data format...")
    Sys.sleep(1)
    setProgress(message=load_files()$mms)
    Sys.sleep(1)
    if (load_files()$mms=="Data format valid!"){
      if (input$Mode=="Fixed tolerance window"){
      setProgress(message="Alignment started...")}
      else{
      setProgress(message="Tuning started...")}
      Sys.sleep(1)
      setProgress(message=alignment()$mms)
    }
  })
  if (load_files()$mms=="Data format valid!"){
    output$blank_message1<-renderText({alignment()$mms})
    
    if (input$enhance){
    output$blank_message2<-renderText({paste0("Enhanced alignment has reduced ",alignment()$loss," molecular features")})
    }
    else{
      output$blank_message2<-renderText({NULL})
      }
  }      
  else {
    updateActionButton(session, "goButton",label = "Try again!")
    output$blank_message1<-renderText({load_files()$mms})
  }
})


## Refresh interface: 
observeEvent(input$killButton,{shinyjs::js$refresh()})


## Change interface:
output$Controls<-renderUI({
  if (input$Mode=="Fixed tolerance window"){
    tagList(
      h3("Mass tolerance window (ppm):"),
      numericInput('mtol', label='',10,min=0,max=30,step=1),
      h3("RT tolerance window (min):"),
      numericInput('rtol', label='',0.2,min=0,max=1,step=0.05),
      h3("Max. percentage of missing values (zeros) allowed:"),
      numericInput('zrt', label='',80,min=0,max=95,step=5),
      checkboxInput("enhance", label = "Enhanced alignment", value = TRUE)
      )}
  else { # For Tune!
    tagList(
      h3("Mass window between (ppm):"),
      sliderInput("range_mtol", label="", min = 0, max = 30, value = c(0,20)),      
      h3("With an increment of (ppm):"),
      numericInput('seg_mtol', label="",3,min=0.5,max=10,step=0.5),
      h3("Retention time window between (min):"),
      sliderInput("range_rtol", label="", min = 0, max = 1, value = c(0,0.5)),      
      h3("With an increment of (min):"),
      numericInput('seg_rtol', label="",0.1,min=0.01,max=0.3,step=0.01))}
})

## Tuning results output
output$heat <- renderD3heatmap({
  x=as.character(alignment()$mtol_test)
  x[length(x)]=paste0(x[length(x)]," ppm")
  y=as.character(alignment()$rtol_test)
  y[length(y)]=paste0(y[length(y)]," min")
  z=alignment()$tune_nbf
  if (!is.null(z)){
#  write.table(z,file="tmp.txt",col.names=F,row.names=F)
    z=data.frame(z)
    colnames(z)=y
    rownames(z)=x
    d3heatmap(x = -z,
            Colv = NULL,
            scale= "none",
            dendrogram="none",
            key = FALSE,
            yaxis_font_size = "12pt",
            xaxis_font_size = "9pt")}
})

output$Tune_choose<-renderUI({
  z=alignment()$tune_nbf
  if (!is.null(z)){
  tagList(
    h3("Choose mass window (ppm):"),
    sliderInput("tune_mtol", label="", min =input$range_mtol[1], max = input$range_mtol[2], value = input$range_mtol[1], step=input$seg_mtol),      
    h3("Which RT window (min):"),
    sliderInput("tune_rtol", label="", min =input$range_rtol[1], max = input$range_rtol[2], value = input$range_rtol[1], step=input$seg_rtol)   
  )}
})

output$Tune_plot<-renderPlot({
  par(mfrow=c(2,2))
  z=alignment()$tune_nbf
  if (!is.null(z)){
    final_list=alignment()$tune_list
    label=paste0(input$tune_mtol,"-",input$tune_rtol)
    final= final_list[[label]]
   if (class(final)!="try-error" & !is.null(final)){
    x=c(final$data[,1],final$data[,1])
    y=c(final$mdev[,1],final$mdev[,2])
    plot(x,y,pch='',ylab='Mass deviation (ppm)',xlab='Mass')
    text(x,y,c(final$mdev[,3],final$mdev[,4]),cex=1.2)
      
    x=c(final$data[,2],final$data[,2])
    y=c(final$rdev[,1],final$rdev[,2])
    plot(x,y,pch='',ylab='Retention time deviation (min)',xlab='Retention time (min)')
    text(x,y,c(final$rdev[,3],final$rdev[,4]),cex=1.2)
    
    results=final$data
    zlist=c()
    for (i in 1:nrow(results)){
      zeros=sum(results[i,3:ncol(results)]==0)
      zlist=c(zlist,zeros)
    }
    hist(zlist,main="",xlab="Number of zeros in each molecular feature")  
  }}
})

## Change interface:
output$Control_fix<-renderUI({
 if (input$Fixed){
  t1=input$Calibrant
  t2=input$End
  tagList(
    sliderInput("range_tt", label="", min = t1, max = t2, value = c(1,3)),      
    h3("Fixed tolerance window (min):"),
    numericInput('tol2', label='',0.2,min=0,max=1,step=0.05)
   )
  }
})


## Visualization tabpanel 

output$Mdev<-renderPlot({ 
  final=alignment()$final
  if (class(final)!="try-error" & !is.null(final)){
    x=c(final$data[,1],final$data[,1])
    y=c(final$mdev[,1],final$mdev[,2])
    plot(x,y,pch='',ylab='Mass deviation (ppm)',xlab='Mass')
    text(x,y,c(final$mdev[,3],final$mdev[,4]),cex=1.2)
  }
}) 

# RT deviation:
output$Rdev<-renderPlot({ 
  final=alignment()$final
  if (class(final)!="try-error" & !is.null(final)){
    x=c(final$data[,2],final$data[,2])
    y=c(final$rdev[,1],final$rdev[,2])
    plot(x,y,pch='',ylab='Retention time deviation (min)',xlab='Retention time (min)')
    text(x,y,c(final$rdev[,3],final$rdev[,4]),cex=1.2)
  }
}) 

# Missing value:
output$Zeros<-renderPlot({ 
  final=alignment()$final
  mini=round(load_files()$NBS*input$zrt/100)
  if (class(final)!="try-error" & !is.null(final)){
    results=final$data
    zlist=c()
    for (i in 1:nrow(results)){
      zeros=sum(results[i,3:ncol(results)]==0)
      zlist=c(zlist,zeros)
    }
  valid=which(zlist<=mini)
  
  h=hist(zlist,plot=FALSE)  
  colors=rep("black",length(h$breaks))
  colors[h$breaks<=mini]="red"
  h=hist(zlist,main="",xlab="Number of zeros in each molecular feature",col=colors)  
  multiTitle(color="black",paste0(nrow(results),"->"),color="red",as.character(length(valid)))
  
  }
})

# ID-Ref:
output$table1 <- renderDataTable({
  filenames=load_files()$filenames
  output=cbind(1:length(filenames),filenames)
  colnames(output)=c("ID","File name")
  return(datatable(output,escape=c(TRUE,TRUE), rownames = F))
})

# Download:

output$downloadMatrix<-downloadHandler(
  filename = function() {
    "Aligned_filtered.txt"},
  content = function(file) {
    final=alignment()$final
    mini=round(load_files()$NBS*input$zrt/100)
    if (class(final)!="try-error" & !is.null(final)){
      results=final$data
      results=results[order(results[,1]),]
      RT_rounded=round(results[,2],1)
      results=results[order(RT_rounded),]
      zlist=c()
      for (i in 1:nrow(results)){
        zeros=sum(results[i,3:ncol(results)]==0)
        zlist=c(zlist,zeros)}
      valid=which(zlist<=mini)
      results=results[valid,]
      write.table(results,file,sep="\t",dec=",",row.names=F,col.names=T)}
  })

## Individual sample evaluation:

# Change start:
 output$Sample_start1<-renderUI({
    NBS=load_files()$NBS
    sliderInput("slider1", label = h3("Starting sample:"), min = 1, max = NBS, value = 1,step=5)
 })

 output$Sample_start2<-renderUI({
   NBS=load_files()$NBS
   sliderInput("slider2", label = h3("Starting sample:"), min = 1, max = NBS, value = 1,step=5)
 })
 
 output$plot_mass_sample<-renderPlot({
   par(mfrow=c(1,5))
   NBS=load_files()$NBS
   filenames=load_files()$filenames
   final=alignment()$final
   if (!is.null(input$slider1)){
   begin=input$slider1
   end=input$slider1+4
   if (end>NBS){end=NBS}
   for (i in begin:end){
      y=final$data[,"mass"]
      x=final$mass_dev_matrix[,i]*1000
      plot(x,y,main=filenames[i],xlab="Deviation (mDa)",ylab="Mass (Da)",xlim=c(-5,5),ylim=c(100,1200))
   }}
 })

 output$plot_RT_sample<-renderPlot({
   par(mfrow=c(1,5))
   NBS=load_files()$NBS
   filenames=load_files()$filenames
   final=alignment()$final
   if (!is.null(input$slider2)){
   begin=input$slider2
   end=input$slider2+4
   if (end>NBS){end=NBS}
   for (i in begin:end){
     y=final$data[,"RT"]
     x=final$RT_dev_matrix[,i]*60
     plot(x,y,main=filenames[i],xlab="Deviation (s)",ylab="Retention time (min)",xlim=c(-15,15),ylim=c(input$Calibrant,input$End))
   }}
 })
 

})
