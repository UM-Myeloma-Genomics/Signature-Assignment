datasetname <- "tMN"  
nsigs <- list(tMN = 3) #number of de novo extracted signatures  
signatures <- list()

cosmDBS<- read.table("COSMIC_v3.2_DBS_GRCh37.txt", sep = "\t", stringsAsFactors = F, header = T) #load cosmic signatures
#add any non-COSMIC signatures desired

SPdbs<- read.table("DBS78_S3_Signatures.txt", sep = "\t", stringsAsFactors = F, header = T) #load sigprofiler de novo extraction
colnames(SPdbs)[1]<- "Type"

#format signatures
for (d in datasetname){ 
  signatures[[d]] <- data.frame(row.names = cosmDBS$Type)   
  for (i in 1:nsigs[[d]]){
    tmpsig <- SPdbs[,-1][i]
    colnames(tmpsig) <- paste0("process_",i)
    signatures[[d]] <- cbind(signatures[[d]],tmpsig) 
  }
}
#pairwise signature assignment
sims<- list()
for (d in datasetname){
  sims[[d]] <- list()
  for (i in 1:nsigs[[d]]){
    for (x in 1:(ncol(cosmDBS[,-1])-1)) {
      for (y in (x+1):(ncol(cosmDBS[,-1]))) {
        a <- cosmDBS[,-1][,c(x,y)]
        r <- nnls(A = as.matrix(a),b = as.vector(signatures[[d]][,i]))
        coeffnorm <- r$x/sum(r$x)
        sims[[d]][[paste0("process_",i)]] <- rbind(sims[[d]][[paste0("process_",i)]],
                                                   data.frame(list(sig = paste0(colnames(a)[1],"+", colnames(a)[2])),   
                                                              cossim=cos.sim2(as.vector(as.matrix(a) %*% r$x),as.vector(signatures[[d]][,i])),
                                                              proportions=paste0("C",x,"=",sprintf("%.2f",coeffnorm[1]),";C",y,"=",sprintf("%.2f",coeffnorm[2])),stringsAsFactors = FALSE))
      }
    }
    sims[[d]][[paste0("process_",i)]] <- sims[[d]][[paste0("process_",i)]][order(sims[[d]][[paste0("process_",i)]]$cossim,decreasing = TRUE),]
  }}

#save tables to file
for (d in datasetname){
  for (i in 1:nsigs[[d]]){
    write.table(sims[[d]][[paste0("process_",i)]],file = paste0("~/cosig.cosine/","DBS", d, "_process_",i,".tsv"),
                sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }
}