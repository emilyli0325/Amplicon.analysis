### Miseq Amplicon-sequencing analysis

# All analysis was performed in mendel. I will use the 12.2019 Shalini run as example to show how to perform this kind of analysis. 

# Three tools will be used in this analysis: bcl2fastq, extract.pl and FASTX-Toolkit.
# bcl2fastq is from illumina for decoding sequencing files.
# extract.pl is a home made perl script. 
# FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) is from Hannon lab. 

############################################################
############ (1) Demultiplex reads by samples  #############
############################################################

# Put the Miseq data in a folder. 
mkdir Analysis
cd /data/infectious/malaria_XUE/Amplicon.CRISPR.Shalini/191218_M01370_0002_000000000-CR6HD

# allow you to open 4000 files at one time
ulimit -n 4000 
# demultiplex files will be output into 'Outs' folder. Decoding summary will be in Outs/Reports folder. Takes about 20 min. 
bcl2fastq -R ./ -o Outs/ 

# move Outs folder to Analysis folder
mv Outs /data/infectious/malaria_XUE/Amplicon.CRISPR.Shalini/Analysis/

# uncompress read files
gunzip * 


############################################################
############ (2) Reverse and complement read2 ##############
############################################################

# You can also use only one read if the location of target sequences is known. Here I used both reads only because I didn't look into the target locations. 

for file in $(ls *R2_001.fastq)
do 
fastx_reverse_complement -i $file -o $file-reco -Q 33
done

## merge two reads
mkdir mergeR1R2

files=$(ls *_L001_R2_001.fastq|sed 's/_L001_R2_001.fastq$//')

for file in ${files[*]}
do
cat $file\_L001_R1_001.fastq $file\_L001_R2_001.fastq-reco > mergeR1R2/$file-merge.fastq
done

# clean up files, it's a personal hobby, I like my folders to be simple and clean. 
mkdir raw
mv *.fastq raw/
mv *-reco raw/

############################################################
############ (3) Count reads carrying different SNPs #######
############################################################

# We use extract.pl to seperate reads carrying different SNPs. It finds the signature sequences and extract them to a specific folder. The signature of target sequences we used here were as shown in Figure 1C (Nair et al., 2018, AAC). You can also define your own sequences depends on different experiment designs.

# use 'perl extract.pl -h' to look at the command options. 
# You can copy the perl script 'extract.pl' to your own folder or use it directly by adding the pathway in front as shown below:

 perl /master/xli/software/HomeMade/extract.pl -h
    i    :input file
    o    :outputfile
    b    :target site : ATCATC[ATGC]{58}TGTTGC
    l    :length of read
    Q1   :threshold for low quality score [20].
    Q2   :maximum percent of low-quality bases[0.2].
    h    :display the help information.

# copy to current folder: cp /master/xli/software/HomeMade/extract.pl ./

## seperate reads. About 1 hour. 
mkdir Anchored Ref 580Mut 580Wt

for file in $(ls mergeR1R2/)
do

perl extract.pl -i mergeR1R2/$file -o Anchored/$file -b ATCATC[ATGC]{58}TGTTGC -l 70
perl extract.pl -i mergeR1R2/$file -o Ref/$file -b ATCATCG[ATGC]{6}A[ATGC]{2}A[ATGC]{26}C[ATGC]{2}T[ATGC]{16}GTGTTGC -l 70
perl extract.pl -i mergeR1R2/$file -o 580Mut/$file -b ATCATCG[ATGC]{6}A[ATGC]{2}A[ATGC]{26}G[ATGC]{2}A[ATGC]{16}ATGTTGC -l 70
perl extract.pl -i mergeR1R2/$file -o 580Wt/$file -b ATCATCG[ATGC]{6}A[ATGC]{2}A[ATGC]{26}G[ATGC]{2}A[ATGC]{16}GTGTTGC -l 70

done

## count reads
for file in $(ls mergeR1R2/)
do

awk '{s++}END{print s/4}' mergeR1R2/$file >> Total.txt
awk '{s++}END{print s/2}' Anchored/$file >> Anchored.txt
awk '{s++}END{print s/2}' 580Mut/$file >> 580Mut.txt
awk '{s++}END{print s/2}' 580Wt/$file >> 580Wt.txt

done

# sample IDs are listed in SampleID.txt.
ls mergeR1R2 > SampleID.txt

## copy numbers in Total.txt,Anchored.txt, 580Mut.txt, 580Wt.txt into an excel(ReadsCount.csv), and add sample information. Remove samples with less than 100 reads covered.

############################################################
############ (4) Plot results with R #######################
############################################################

# Prepare input files

setwd("C:/Users/xli/Documents/Emily/CRISPR/3_Amplicon_12.2019/Output")

rawdata2 <- read.csv("ReadsCount.C580Y.clean.csv", header=TRUE, sep=",", check.names=FALSE, row.names=1)
A.1 <-rawdata2[grep("N", rawdata2$CultureID),]
B.1 <-rawdata2[grep("M", rawdata2$CultureID),]
C.1 <-rawdata2[grep("S", rawdata2$CultureID),]


# Figure 2 selection coefficient

pdf('Figure2_C580Y_clean.pdf', width=15, height=10)
layout(rbind(c(1,2,3),c(4,5,6)))
par(oma = c(0.5,0.5,0.5,0.5), mai= c(0.6,0.6,0.5,0.2))


# A
plot(A.1$Asexual.cycle, A.1$Percent.580Mut, main = "N", sub="", xlab="48h asexual cycles",  ylab = "580MT proportion", cex.lab = 2, cex.axis = 1.8, cex.main = 2, ylim=range(0:1),pch=21, col= "black", bg = c("red", "blue", "green", "black", "Magenta", "Cyan")[as.numeric(A.1$RepeatID)], cex=2)
# B
plot(B.1$Asexual.cycle, B.1$Percent.580Mut, main = "M", sub="", xlab="48h asexual cycles",  ylab = "580MT proportion", cex.lab = 2, cex.axis = 1.8, cex.main = 2, ylim=range(0:1),pch=21, col= "black", bg = c("red", "blue", "green", "black", "Magenta", "Cyan")[as.numeric(A.1$RepeatID)], cex=2)
# C
plot(C.1$Asexual.cycle, C.1$Percent.580Mut, main = "S", sub="", xlab="48h asexual cycles",  ylab = "580MT proportion", cex.lab = 2, cex.axis = 1.8, cex.main = 2, ylim=range(0:1),pch=21, col= "black", bg = c("red", "blue", "green", "black", "Magenta", "Cyan")[as.numeric(A.1$RepeatID)], cex=2)

# D
plot(A.1$Asexual.cycle, log(A.1$Percent.580Mut/A.1$Percent.580Wt), main = "N", sub="", xlab="48h asexual cycles",  ylab = "ln(genotype ratio)", ylim=range(-8:1), cex.lab = 2, cex.axis = 1.8, cex.main = 2, pch=21, col= c("red", "blue", "green", "black", "Magenta", "Cyan")[as.numeric(A.1$RepeatID)], bg = NA, cex=2)
A.mat <- cbind(A.1[,c(6,7)],log(A.1$Percent.580Mut/A.1$Percent.580Wt))
colnames(A.mat)[3] <- "LN"
A.mat <- A.mat[is.finite(rowSums(A.mat)),]
m.1 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "1")
m.2 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "2")
m.3 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "3")
m.4 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "4")
m.5 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "5")
m.6 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "6")
abline(m.1, col = "red", cex=2)
abline(m.2, col = "blue", cex=2)
abline(m.3, col = "green", cex=2)
abline(m.4, col = "black", cex=2)
abline(m.5, col = "Magenta", cex=2)
abline(m.6, col = "Cyan", cex=2)

# E
plot(B.1$Asexual.cycle, log(B.1$Percent.580Mut/B.1$Percent.580Wt), main = "M", sub="", xlab="48h asexual cycles",  ylab = "ln(genotype ratio)", ylim=range(-8:1),cex.lab = 2, cex.axis = 1.8, cex.main = 2, pch=21, col= c("red", "blue", "green", "black", "Magenta", "Cyan")[as.numeric(B.1$RepeatID)], bg = NA, cex=2)
B.mat <- cbind(B.1[,c(6,7)],log(B.1$Percent.580Mut/B.1$Percent.580Wt))
colnames(B.mat)[3] <- "LN"
B.mat <- B.mat[is.finite(rowSums(B.mat)),]
m.1 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "1")
m.2 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "2")
m.3 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "3")
m.4 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "4")
m.5 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "5")
m.6 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "6")
abline(m.1, col = "red", cex=2)
abline(m.2, col = "blue", cex=2)
abline(m.3, col = "green", cex=2)
abline(m.4, col = "black", cex=2)
abline(m.5, col = "Magenta", cex=2)
abline(m.6, col = "Cyan", cex=2)


# F
plot(C.1$Asexual.cycle, log(C.1$Percent.580Mut/C.1$Percent.580Wt), main = "S", sub="", xlab="48h asexual cycles",  ylab = "ln(genotype ratio)", ylim=range(-8:1), cex.lab = 2, cex.axis = 1.8, cex.main = 2, pch=21, col= c("red", "blue", "green", "black", "Magenta", "Cyan")[as.numeric(C.1$RepeatID)], bg = NA, cex=2)
C.mat <- cbind(C.1[,c(6,7)],log(C.1$Percent.580Mut/C.1$Percent.580Wt))
colnames(C.mat)[3] <- "LN"
C.mat <- C.mat[is.finite(rowSums(C.mat)),]
m.1 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "1")
m.2 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "2")
m.3 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "3")
m.4 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "4")
m.5 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "5")
m.6 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "6")
abline(m.1, col = "red", cex=2)
abline(m.2, col = "blue", cex=2)
abline(m.3, col = "green", cex=2)
abline(m.4, col = "black", cex=2)
abline(m.5, col = "Magenta", cex=2)
abline(m.6, col = "Cyan", cex=2)

dev.off()


# Table 1 significant test. Most of these are repeating Figure 2, I put them here for better understanding. 
# put AssayID,s,se,R2,n,P IN SLOPE.csv


A.mat <- cbind(A.1[,c(6,7)],log(A.1$Percent.580Mut/A.1$Percent.580Wt))
colnames(A.mat)[3] <- "LN"
A.mat <- A.mat[is.finite(rowSums(A.mat)),]
m.1 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "1")
m.2 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "2")
m.3 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "3")
m.4 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "4")
m.5 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "5")
m.6 <- lm(LN~Asexual.cycle, data = A.mat, subset = RepeatID == "6")
summary(m.1)
summary(m.2)
summary(m.3)
summary(m.4)
summary(m.5)
summary(m.6)
confint(m.1, level=0.95)
confint(m.2, level=0.95)
confint(m.3, level=0.95)
confint(m.4, level=0.95)
confint(m.5, level=0.95)
confint(m.6, level=0.95)

m <- lm(LN~Asexual.cycle, data = A.mat)
summary(m)
confint(m, level=0.95)


B.mat <- cbind(B.1[,c(6,7)],log(B.1$Percent.580Mut/B.1$Percent.580Wt))
colnames(B.mat)[3] <- "LN"
B.mat <- B.mat[is.finite(rowSums(B.mat)),]
m.1 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "1")
m.2 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "2")
m.3 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "3")
m.4 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "4")
m.5 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "5")
m.6 <- lm(LN~Asexual.cycle, data = B.mat, subset = RepeatID == "6")
summary(m.1)
summary(m.2)
summary(m.3)
summary(m.4)
summary(m.5)
summary(m.6)
confint(m.1, level=0.95)
confint(m.2, level=0.95)
confint(m.3, level=0.95)
confint(m.4, level=0.95)
confint(m.5, level=0.95)
confint(m.6, level=0.95)

m <- lm(LN~Asexual.cycle, data = B.mat)
summary(m)
confint(m, level=0.95)



C.mat <- cbind(C.1[,c(6,7)],log(C.1$Percent.580Mut/C.1$Percent.580Wt))
colnames(C.mat)[3] <- "LN"
C.mat <- C.mat[is.finite(rowSums(C.mat)),]
m.1 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "1")
m.2 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "2")
m.3 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "3")
m.4 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "4")
m.5 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "5")
m.6 <- lm(LN~Asexual.cycle, data = C.mat, subset = RepeatID == "6")
summary(m.1)
summary(m.2)
summary(m.3)
summary(m.4)
summary(m.5)
summary(m.6)
confint(m.1, level=0.95)
confint(m.2, level=0.95)
confint(m.3, level=0.95)
confint(m.4, level=0.95)
confint(m.5, level=0.95)
confint(m.6, level=0.95)

m <- lm(LN~Asexual.cycle, data = C.mat)
summary(m)
confint(m, level=0.95)



# Figure 3 summary of head-to-head competition

# merge result from different replicates 
library(metafor)
library(metap)

data <- read.csv("SLOPES.csv", header=TRUE, sep=",", check.names=FALSE)
Comp1 <-data[grep("N", data$CompetitionID),]
Comp2 <-data[grep("M", data$CompetitionID),]
Comp3 <-data[grep("S", data$CompetitionID),]
meta1<-rma(s, se, data=Comp1, measure="OR")
meta2<-rma(s, se, data=Comp2, measure="OR")
meta3<-rma(s, se, data=Comp3, measure="OR")

meta1
meta2
meta3

sumlog(Comp1$P)
sumlog(Comp2$P)
sumlog(Comp3$P)

# forest figure
library(ggplot2)

data <- read.csv("forest.csv", header=TRUE, sep=",", check.names=FALSE)
# reverses the factor level ordering for labels after coord_flip()
data$label <- factor(data$label, levels=rev(data$label))
pdf('Figure3.pdf', width=4, height=6)
par(oma = c(0.5,0.5,0.5,0.5), mai= c(0.6,0.6,0.5,0.2))
# MAKE INITIAL PLOT
fp <- ggplot(data=data, aes(x=label, y=s, ymin=lower.CI, ymax=upper.CI)) +
        geom_pointrange(shape=21, size=0.8, col= "black", bg=c("grey", "red")[as.numeric(data$col)]) + 
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
		scale_x_discrete(labels = rev(data$AssayID)) +
        theme_bw()  # use a white background
# ADD AXIS LABELS3
fp <- fp + ylab("Selection coefficient (95% CI)") + xlab("Assay ID") + ggtitle("") + theme(legend.position="none")
print(fp)
dev.off()
