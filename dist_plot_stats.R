#library(readr)

args<-commandArgs(TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied (sim input file), (real input file), (header)", call.=FALSE)
} else{
  print(args)
}

file_sim<-args[1] #simulation output, should have the form input_ABCtoolbox_M1_HPC.txt
file_real<-args[2] #real data, should have the form real_output_M23.summary
header<-args[3] #header of simulation containing desired columns 

input_ABCtoolbox<-read.table(file_sim, header = T);
real_output<-read.table(file_real, header = T);
keep<-scan(header, character(), quote="")

out_params<-paste(strsplit(file_sim, ".txt")[[1]],"_params.pdf", sep="")
out_stats<-paste(strsplit(file_sim, ".txt")[[1]],"_stats.pdf", sep="")
out_pca<-paste(strsplit(file_sim, ".txt")[[1]],"_pca.pdf", sep="")
print(out_params)
print(out_stats)
print(out_pca)

#input_ABCtoolbox<-read_delim("/vol_c/results_macsSwig_AJmodels_instant/input_ABCtoolbox_M1_HPC.txt","\t", escape_double = FALSE, trim_ws = TRUE)
#real_output <- read_delim("/vol_c/ABC_AJmodels/real_output_23.summary","\t", escape_double = FALSE, trim_ws = TRUE)
#keep<-scan("/vol_c/results_macsSwig_AJmodels_instant/header_M1_132.txt", character(), quote="")

input<-input_ABCtoolbox[keep]
n_params<-which(colnames(input)=="SegS_Af_CGI")

#names(input[,1:n_params]) <- c("sim","Asc_NAF","Asc_NEU","Asc_NCHB","daf","Log10_NAF","Log10_NANC","Log10_NCEU","Log10_NCHB","Log10_NWA","Log10_NEA","Log10_NJ","Log10_NM","m","Tgrowth_Af","TAF","TEM","Teu_as","TA","TMJ","TAEW","Tm","rWA","rEA","rMJ")

ss=names(input)
Tbl <- as.data.frame(input)
obs<-Tbl[1:nrow(Tbl),n_params:ncol(Tbl)]
removals<-c()
for(i in 1:ncol(obs)){
  if(all(obs[,i] == obs[[1,i]])){
    removals<-append(removals, i)
  }
}
if(length(removals)==0){
  new.obs <- obs
} else {
  new.obs <- obs[,-removals] 
}
rm(obs)

## Plot distribution of parameters
pdf(out_params, width=12, height=9);
par(mfrow=c(3,4));
for (j in 2:(n_params-1)) {
  plot(density(Tbl[ , j]), main=names(Tbl)[j])
}
dev.off();

## Plot distribution of summary statistics
pdf(out_stats, width=12, height=9);
par(mfrow=c(3,4));
for (j in 1:ncol(new.obs)) {
  plot(density(new.obs[ , j]), main=names(new.obs)[j])
  abline(v = real_output[,j][1], col = "red")
}
dev.off();


#PCA of the simulations

##calculate the PCs based just on the simulations
mydata.pca <- prcomp(new.obs, retx=TRUE, center=TRUE,scale.=TRUE)
scores <-mydata.pca$x
sd <- mydata.pca$sdev
eigenvalues <- sd^2

##project the observed data into the PC of the simulations
projection<-scale(real_output,mydata.pca$center,mydata.pca$scale) %*% mydata.pca$rotation

##find min and max to make the plots
min1<-min(scores[,1])
max1<-max(scores[,1])
min2<-min(scores[,2])
max2<-max(scores[,2])
if (projection[1] < min1){ min1<-projection[1] }
if (projection[1] > max1){ max1<-projection[1] }
if (projection[2] < min2){ min2<-projection[2] }
if (projection[2] > max2){ max2<-projection[2] }

###plot
pdf(out_pca, width=6, height=6);
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2],pch=20,cex=0.2,col="white", xlim=c(min1,max1), ylim=c(min2,max2), xlab=xlab, ylab=ylab)
points(scores[,1][2:1001], scores[,2][2:1001],col="grey",pch=20)
points(projection[1],projection[2],col="red",pch=20)
dev.off();