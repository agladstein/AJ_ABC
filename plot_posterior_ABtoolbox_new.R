##Script to plot the posterior distributions, refactored from a script by Consuelo Quinto, which was modified from a script provided in the first version of ABCtoolbox
##Author: Ariella Gladstein
##August 17, 2017

##Rscript plot_posterior_ABtoolbox_new.R test results_regions_5percent_ss_10kb_Affy6.0_snps_transformer_20pls.txt_100_model0_MarginalPosteriorDensities_Obs0.txt results_regions_5percent_ss_10kb_Affy6.0_snps_transformer_20pls.txt_100_model0_BestSimsParamStats_Obs0.txt hpdi_values; mv out.pdf plot_posterior_densities_20pls_100retsims.pdf

## Rscript /vol_c/src/macsswig_simsaj/plot_posterior_ABtoolbox_new.R keepPowerStats_input_ABCtoolbox_M2_HPC_OSG_2.txt ABC_correlatedstats6_1446125_pruneCorStats_90_model0_MarginalPosteriorDensities_Obs0.txt ABC_correlatedstats6_1446125_pruneCorStats_90_model0_BestSimsParamStats_Obs0.txt

####################
args<-commandArgs(TRUE)
if (length(args)<4) {
  stop("Four arguments must be supplied (parameter values from all sims), (real input file), (header), (out file)", call.=FALSE)
} else{
  print(args)
}

simfile<-args[1] ##parameter values from all sims. Simulation output, should have the form input_ABCtoolbox_M1_HPC.txt
posfile<-args[2] ##posterior distribution, *MarginalPosteriorDensities_Obs0.txt
retfile<-args[3] ##retained simulations *BestSimsParamStats_Obs0.txt
#hpdivals<-args[4] ##file with the 95 HPDI values

sims<-read.table(simfile, header=T, nrows=20000);
post<-read.table(posfile, header=T);
ret<-read.table(retfile, header=T);
#hpdi<-read.table(hpdivals,header=T, row.names=1);

out<-args[4]
print(out)

names_params<-c("Asc_NAF","Asc_NEU","Asc_NCHB","daf", "Log10_NAF","Log10_NCEU","Log10_NCHB","Log10_NWA","Log10_NEA","Log10_NAg","Log10_NJ","Log10_NM","m","TAF","TEM","Teu_as","TA","TMJ","TAEW","Tm");
params<-c(2,3,4,5,6,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23);

pdf(out, width=9, height=12);
par(mfrow=c(5,4));

for(j in 1:length(params)){
  
  p<-density(sims[,params[j]]);
  xmax<-max(p$x);
  #print(xmax);
  xmin<-min(p$x);
  ymax<-max(p$y);
  ymax<-max(max(p$y), max(post[,2*j+1]))*1.1;
  #print(params[j])
  #print(hpdi[1,j])
  #print(hpdi[params[j]])
  
  plot(post[,2*j], post[,2*j+1], main=names_params[j], xlab=paste("mode at ", post[post[,2*j+1]==max(post[,2*j+1]),2*j]), ylab="Density", type='l', col="red", xlim=c(xmin, xmax), ylim=c(0,ymax));
  lines(p, col="black");
  lines(density(ret[,j+2]), col="blue");
  #abline(v = hpdi[1,j], col = "dimgrey", lwd=2.5, lty=2);
  #abline(v = hpdi[2,j], col = "dimgrey", lwd=2.5, lty=2);
}
dev.off();
