#!/bin/sh

#NOTE: run this script from within the directory containing the QC'ed files
#NOTE: this file must be run after QC.sh as it uses the trimmed_pruned files
echo "Have you run QC.sh? Respond Y or N, followed by [ENTER]: "
read QCdone
#NOTE: there must be an individual list with FID and IID columns (no column labels needed)
echo "Do you have an individual list with two columns (FID and IID) in this folder? Respond Y or N, followed by [ENTER]: "
read indlist

if [ $QCdone == "Y" ] && [ $indlist == "Y" ]
then
	echo "good to go"
else
	echo "You need to prepare the files for this script to run properly"
fi

module load plink
module load r

echo "Enter name of study population (e.g. WSC, MrOS, APOE), followed by [ENTER]: "
read study

echo "Enter number of genotype/chip/array pseudocohorts, followed by [ENTER]: "
read cohortnum

# gather genotype/chip/array pseudocohort(s) names into an array called "list"
for i in {1..$cohortnum}
do
	echo "Enter name(s) of genotype/chip/array pseudocohort(s), separated by spaces, followed by [ENTER]: "
	read cohortnames
	list=($cohortnames)
done

# determine phenotype name (as designated in the ${cohortname}_pheno.txt file)
#	can generate the pheno file, using the fam2pheno.R script: https://www.dropbox.com/s/g8e5zzwkvdc8ny0/fam2pheno.R?dl=0
echo "Enter BINARY phenotype name as it appears in the {cohortname}_pheno.txt file, followed by [ENTER]: "
read pheno

workDir=$PWD

# loop over cohort names in "list" to submit jobs
for k in "${list[@]}"
do

	# make directory to receive files for pseudocohort output
	mkdir -p $workDir/${k}_subset
	cp ${k}_trimmed_pruned.* ${k}_subset
	#generate study-of-interest only subset
	plink --file ${k}_trimmed_pruned --keep ${study}.indlist --allow-no-sex --make-bed --out $workDir/${k}_subset/${k}_${study}only
	#get SNPlists of each of the cohort's PLINK BINARIES
	plink --bfile $workDir/${k}_subset/${k}_${study}only --allow-no-sex --write-snplist --make-bed --out $workDir/${k}_subset/${k}
	#do a "cross extraction" of SNPs from HapMap and the cohorts
	plink --bfile /srv/gsfs0/projects/mignot/PLMGWAS/HapMap4QC/hapmap3_r2_b37_fwd.consensus.qc.poly --allow-no-sex --extract $workDir/${k}_subset/${k}.snplist --make-bed --out $workDir/${k}_subset/hapmap_${k}extract
	plink --bfile $workDir/${k}_subset/${k} --allow-no-sex --extract /srv/gsfs0/projects/mignot/PLMGWAS/HapMap4QC/hapmap3_r2_b37_fwd.consensus.qc.poly.snplist --make-bed --out $workDir/${k}_subset/${k}_HPextract
	#join the SNP lists
	plink --bfile $workDir/${k}_subset/${k}_HPextract --bmerge $workDir/${k}_subset/hapmap_${k}extract.bed $workDir/${k}_subset/hapmap_${k}extract.bim $workDir/${k}_subset/hapmap_${k}extract.fam --allow-no-sex --make-bed --out $workDir/${k}_subset/${k}_HPmerge
	count=`ls -1 ${k}_subset/*.missnp 2>/dev/null | wc -l`
	if [ $count != 0 ]
	then 
	plink --bfile $workDir/${k}_subset/hapmap_${k}extract --exclude $workDir/${k}_subset/${k}_HPmerge-merge.missnp --make-bed --out $workDir/${k}_subset/hapmap_${k}extract_SNPmis
	plink --bfile $workDir/${k}_subset/${k}_HPextract --exclude $workDir/${k}_subset/${k}_HPmerge-merge.missnp --make-bed --out $workDir/${k}_subset/${k}_HPextract_SNPmis
	plink --bfile $workDir/${k}_subset/${k}_HPextract_SNPmis --bmerge $workDir/${k}_subset/hapmap_${k}extract_SNPmis.bed $workDir/${k}_subset/hapmap_${k}extract_SNPmis.bim $workDir/${k}_subset/hapmap_${k}extract_SNPmis.fam --allow-no-sex --make-bed --out $workDir/${k}_subset/${k}_HPmerge
	fi
	#next will generate multidimensional scaling analysis
	plink --bfile $workDir/${k}_subset/${k}_HPmerge --allow-no-sex --genome --out $workDir/${k}_subset/${k}_HPmerge_genome
	plink --bfile $workDir/${k}_subset/${k}_HPmerge --allow-no-sex --read-genome $workDir/${k}_subset/${k}_HPmerge_genome.genome --cluster --mds-plot 10 --out $workDir/${k}_subset/${k}_HPmerge_MDS
	#format files from variable-spaces to tab-delimited
	cat $workDir/${k}_subset/${k}_HPmerge_MDS.mds | sed 's/ \+/\t/g' | sed -e 's/^[ \t]*//' > $workDir/${k}_subset/${k}_MDS_space2tab.txt
cat <<EOF >$workDir/${k}_subset/${k}_plotprep.R
tabs <- read.table("${workDir}/${k}_subset/${k}_MDS_space2tab.txt", header = T, sep = "\t")
tabs <- tabs[,c(1:13)]
#this creates a cc designation to differentiate study individuals from HapMap individuals
tabs <- cbind(cc = 0,tabs)
indlist <- read.table("${workDir}/${study}.indlist")
classytabs <- tabs
#this generates a designation of 1 if the individual is in the study of interest, and 0 if from HapMap
classytabs[which(tabs[,"IID"] %in% indlist[,"V2"]),"cc"] <- 1
write.table(classytabs,file="${workDir}/${k}_subset/${k}_tablimited.txt",quote=F,row.names=F,sep="\t")
q()
EOF

	R CMD BATCH $workDir/${k}_subset/${k}_plotprep.R
	
	cat <<EOF >$workDir/${k}_subset/${k}_plot.R
require(graphics)

data <- read.table("${workDir}/${k}_subset/${k}_tablimited.txt",header=TRUE,sep='\t')
ethnicities <- read.table('/srv/gsfs0/projects/mignot/PLMGWAS/PlotAgainstHapMap/ethnicities_hapmap.txt',header=TRUE,sep='\t')

data_hapmap <- data[data[,"cc"]==0,]
data_hapmap_eth <- merge(x = data_hapmap, y = ethnicities, by = 'IID', all = TRUE)
data_hapmap_PC1 <- data_hapmap_eth[,"C1"]
data_hapmap_PC2 <- data_hapmap_eth[,"C2"]
data_hapmap_PC3 <- data_hapmap_eth[,"C3"]
data_${pheno} <- data[data[,"cc"]==1,]
data_${pheno}_PC1 <- data_${pheno}[,"C1"]
data_${pheno}_PC2 <- data_${pheno}[,"C2"]
data_${pheno}_PC3 <- data_${pheno}[,"C3"]

#assess outliers
mean1 <- mean(data_${pheno}_PC1)
stdev1 <- sd(data_${pheno}_PC1)
UL1 <- mean1+3*stdev1
LL1 <- mean1-3*stdev1
outlier1 <- data[which(data[,"C1"]>UL1 | data[,"C1"]<LL1),c("cc","FID","IID","C1")]
outliersub1 <- outlier1[outlier1[,"cc"]==1,c("FID","IID")]
write.table(outliersub1,file="$workDir/${k}_subset/${k}_PC1outliers.txt",quote=F,row.names=F)
mean2 <- mean(data_${pheno}_PC2)
stdev2 <- sd(data_${pheno}_PC2)
UL2 <- mean2+3*stdev2
LL2 <- mean2-3*stdev2
outlier2 <- data[which(data[,"C2"]>UL2 | data[,"C2"]<LL2),c("cc","FID","IID","C2")]
outliersub2 <- outlier2[outlier2[,"cc"]==1,c("FID","IID")]
write.table(outliersub2,file="$workDir/${k}_subset/${k}_PC2outliers.txt",quote=F,row.names=F)
mean3 <- mean(data_${pheno}_PC3)
stdev3 <- sd(data_${pheno}_PC3)
UL3 <- mean3+3*stdev3
LL3 <- mean3-3*stdev3
outlier3 <- data[which(data[,"C3"]>UL3 | data[,"C3"]<LL3),c("cc","FID","IID","C3")]
outliersub3 <- outlier3[outlier3[,"cc"]==1,c("FID","IID")]
write.table(outliersub3,file="$workDir/${k}_subset/${k}_PC3outliers.txt",quote=F,row.names=F)

palette()
palette(rainbow(11))
#${k} data
all_PC1 <- append(data_hapmap_PC1,data_${pheno}_PC1)
all_PC2 <- append(data_hapmap_PC2,data_${pheno}_PC2)
all_PC3 <- append(data_hapmap_PC2,data_${pheno}_PC3)

#plot first 3 PCs to look for outliers
options(bitmapType='cairo')
png("$workDir/${k}_subset/${study}_${k}_1vs2.png", width = 8, height = 8, units = 'in', bg = "transparent", res = 300)
plot(data_${pheno}_PC1,data_${pheno}_PC2,col="white",xlab="PC1",ylab="PC2",main="${k} PC1 vs PC2")
points(data_${pheno}_PC1,data_${pheno}_PC2,col=rgb(0,0,0,1),pch=20,xlab="PC1",ylab="PC2",main="${pheno} PC1 vs PC2")
dev.off()
png("$workDir/${k}_subset/${study}_${k}_1vs3.png", width = 8, height = 8, units = 'in', bg = "transparent", res = 300)
plot(data_${pheno}_PC1,data_${pheno}_PC3,col="white",xlab="PC1",ylab="PC3",main="${k} PC1 vs PC3")
points(data_${pheno}_PC1,data_${pheno}_PC3,col=rgb(0,0,0,1),pch=20,xlab="PC1",ylab="PC3",main="${pheno} PC1 vs PC3")
dev.off()
png("$workDir/${k}_subset/${study}_${k}_2vs3.png", width = 8, height = 8, units = 'in', bg = "transparent", res = 300)
plot(data_${pheno}_PC2,data_${pheno}_PC3,col="white",xlab="PC2",ylab="PC3",main="${k} PC2 vs PC3")
points(data_${pheno}_PC2,data_${pheno}_PC3,col=rgb(0,0,0,1),pch=20,xlab="PC2",ylab="PC3",main="${pheno} PC2 vs PC3")
dev.off()

#plot against HapMap
png("$workDir/${k}_subset/${study}_${k}_against_HapMap.png", width = 8, height = 8, units = 'in', bg = "transparent", res = 300)
plot(all_PC1,all_PC2,col="white",xlab="PC1",ylab="PC2",main="${pheno} PC1 vs PC2")
points(data_hapmap_PC1,data_hapmap_PC2,xlab="PC1",ylab="PC2",main="${pheno} PC1 vs PC2", col=data_hapmap_eth[,"population"])
points(data_${pheno}_PC1,data_${pheno}_PC2,col=rgb(0,0,0,0.5),pch=20,xlab="PC1",ylab="PC2",main="${pheno} PC1 vs PC2")
legend("topleft",legend=levels(factor(data_hapmap_eth[,"population"])),text.col=seq_along(levels(factor(data_hapmap_eth[,"population"]))))
dev.off()
q()
EOF
	
	R CMD BATCH $workDir/${k}_subset/${k}_plot.R
	#You should check this plot to verify that the samples are of the expected ethnicity (i.e. clustered where you expect them)
	#Individuals falling >3 SD beyond the mean in each of the 3 main dimensions are considered outliers
	# -can use --remove in association with {cohortname}_PC*outliers.txt to eliminate if desired
	
	#NOTE: the flag is thrown because of the phenotype assumption: 0=unaffected 1=affected -9=missing
	#Notice the oddity of having to call a specific instance of PLINK, the reasons are below:
	# -at the time of script generation, the SCG instance of PLINK (v1.90b3c 64-bit (2 Feb 2015)) resulted in an error: Zero valid tests; --adjust skipped
	# -this issue resolved upon using a more updated version of PLINK (v1.90b3.32 64-bit (24 Feb 2016)), located in /srv/gsfs0/projects/mignot/PLMGWAS/
	cp /srv/gsfs0/projects/mignot/PLMGWAS/plink1.9 $workDir
	rsync -aq $workDir/${k}_pheno.txt $workDir/${k}_subset/${k}_pheno.txt
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --logistic --adjust --out $workDir/${k}_subset/${k}_unadj
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1 --logistic --adjust --out $workDir/${k}_subset/${k}_C1
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C2 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C2
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C3 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C3
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C4 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C4
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C5 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C5
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C6 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C6
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C7 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C7
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C8 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C8
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C9 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C9
	./plink1.9 --bfile $workDir/${k}_subset/${k}_${study}only --pheno $workDir/${k}_subset/${k}_pheno.txt --pheno-name ${pheno} --covar $workDir/${k}_subset/${k}_HPmerge_MDS.mds --covar-name C1-C10 --logistic --adjust --out $workDir/${k}_subset/${k}_C1-C10
	
	cat <<EOF >$workDir/${k}_subset/${k}_labmdas.txt
Adjustment lambda
EOF
	printf 'unadj ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_unadj.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C2 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C2.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C3 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C3.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C4 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C4.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C5 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C5.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C6 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C6.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C7 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C7.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C8 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C8.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C9 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C9.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	printf 'C1-C10 ' >> $workDir/${k}_subset/${k}_labmdas.txt
	grep 'lambda' ${k}_subset/${k}_C1-C10.log | sed -n 's/.*(based on median chisq) = //p' >> $workDir/${k}_subset/${k}_labmdas.txt
	
	
	cat <<EOF >$workDir/${k}_subset/${k}_MDSselect.R
#determine optimal number of components, based on minimum lambda
lambdas <- read.table("${k}_subset/${k}_labmdas.txt",header=T,sep=" ")
dimnum <- lambdas
dimnum[,"Adjustment"] <- c(0:10)
options(bitmapType='cairo')
png("$workDir/${k}_subset/deflation.png",height=1000,width=1000)
plot(dimnum, xlab="Number of MDS dimensions adjusted for", ylab="lambda", main="Optimal MDS adjustment for ${k}", type="b")
dev.off()
min<-lambdas[which.min(lambdas[,"lambda"]),]
index<-which.min(lambdas[,"lambda"])
#if (index != 1) {
#	lmean<-mean(lambdas[,"lambda])
#	lsd<-sd(lambdas[,"lambda])
#	if (lambdas[index,"lambda"] < lmean+lsd)
#		optimalMDS <- lambdas[which.min(lambdas[,"lambda])-1,]
#		oneless<-TRUE
#	}
#	else {
#	optimalMDS <- lambdas[which.min(lambdas[,"lambda]),]
#	oneless<-FALSE
#}
write.table(min,file="$workDir/${k}_subset/${k}_optimalMDSnumber.txt",quote=F,row.names=F)
chooseMDS <- read.table("${k}_subset/${k}_MDS_space2tab.txt",header=T)
chosenMDS <- chooseMDS[,c(1,2)]
if (index > 1) {
	chosenMDS <- cbind(chosenMDS,chooseMDS[,4:(index+3)])
	write.table(chosenMDS, file="$workDir/${k}_subset/${k}_optimalMDScovs.txt",quote=F,row.names=F)

	#next will generate QQplots
	options(bitmapType='cairo')
	png("$workDir/${k}_subset/qqplot_compare.png", height=1000, width=1000)
	par(mfrow=c(2,1))
	#QQplot unadjusted
	${pheno}_unadj<-read.table("${k}_subset/${k}_unadj.assoc.logistic", header=TRUE)
	${pheno}_unadj.add.p<-${pheno}_unadj[${pheno}_unadj[,"TEST"]==c("ADD"),"P"]
	observed <- sort(${pheno}_unadj.add.p)
	lobs <- -(log10(observed))
	expected <- c(1:length(observed))
	lexp <- -(log10(expected / (length(expected)+1)))
	plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,max(lobs)), ylim=c(0,max(lobs)), las=1, xaxs="i", yaxs="i", bty="l", main = "${pheno} Unadjusted")
	points(lexp, lobs, pch=23, cex=.4, bg="black")

	#QQplot adjusted
	${pheno}_optimalMDS<-read.table(paste("${k}_subset/${k}_",lambdas[index,"Adjustment"],".assoc.logistic",sep=""), header=TRUE)
	${pheno}_optimalMDS.add.p<-${pheno}_optimalMDS[${pheno}_optimalMDS[,"TEST"]==c("ADD"),"P"]
	observed <- sort(${pheno}_unadj.add.p)
	lobs <- -(log10(observed))
	expected <- c(1:length(observed))
	lexp <- -(log10(expected / (length(expected)+1)))
	plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,max(lobs)), ylim=c(0,max(lobs)), las=1, xaxs="i", yaxs="i", bty="l", main = paste("${pheno}",as.character(lambdas[index,"Adjustment"]),sep=" "))
	points(lexp, lobs, pch=23, cex=.4, bg="black")
	dev.off()
} else {
	write.table(chosenMDS, file="$workDir/${k}_subset/${k}_optimalMDScovs.txt",quote=F,row.names=F)
	
	#next will generate QQplots
	options(bitmapType='cairo')
	png("$workDir/${k}_subset/qqplot_no_adjustment_indicated.png", height=1000, width=1000)
	#QQplot unadjusted
	${pheno}_unadj<-read.table("${k}_subset/${k}_unadj.assoc.logistic", header=TRUE)
	${pheno}_unadj.add.p<-${pheno}_unadj[${pheno}_unadj[,"TEST"]==c("ADD"),"P"]
	observed <- sort(${pheno}_unadj.add.p)
	lobs <- -(log10(observed))
	expected <- c(1:length(observed))
	lexp <- -(log10(expected / (length(expected)+1)))
	plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,max(lobs)), ylim=c(0,max(lobs)), las=1, xaxs="i", yaxs="i", bty="l", main = "${pheno} Unadjusted")
	points(lexp, lobs, pch=23, cex=.4, bg="black")
	dev.off()
}

q()
EOF

	R CMD BATCH $workDir/${k}_subset/${k}_MDSselect.R
	
	###WHY DO THE QQplots look so deflated?  Compare to the GenABEL method
	###NEXT STEPS: AND UPDATE the "What's in here.txt" files with what files are useful for next steps, decision points, and publication
	###-make fam2cov.R script to generate the dummy {cohortname}_cov.txt file
	###-fill ${cohortname}_cov.txt with dummy (1/0) variables for pseudocohort, plate number (non dummy)
	###-add in relevant MDS dimensions from optimalMDS.txt
	###-generate script for production of SNPtest .sample files from PLINK cov/pheno files
	
	###-create cut point for partitioning data by ethnicities after looking at MDS plotted against HapMap
	###---can then re-do MDS to account population stratification BY enthic group within each study's pseudocohort/array subset
	###-pipeline into SHAPEIT (pre-phase) and impute2
	###-QC again with EasyQC (another flow diagram and QC.txt file)
	###-PLINK --oxford to get gen/sample files; add cov/pheno to sample files
	###-SNPtest each study subset
	###-META the results
	###-PRSice analysis
	###-VEGAS-STAMS analysis
	
done
