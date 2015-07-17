\name{GBSPipeline}
\alias{GBSPipeline}
\title{
Genotyping-by-Sequencing Pipeline
}
\description{
Uses the 'rJava' package to run Java code for the Genotyping-by-Sequencing pipeline described by Poland et al. (2012). The pipeline makes use of open source Tassel 4.0 software. See \url{http://tassel.bitbucket.org/TasselArchived.html} for how to obtain Tassel 4.0 and documentation about the program.
}
\usage{
GBSPipeline(projectName, keyFile, seqDir, tasselJar, 
heapSize, enzyme = "PstI-MspI", maxMissingData = 0.8, 
minorAlleleFreq = 0.01, maxHet = 0.1, isDHpopulation = F, 
callHets = T, pVal = 0.001, runPipeline = T, isQseq = F, 
tagCountsDir = NULL, masterTagsFile = NULL, tbtDir = NULL, 
tbtFile = NULL, tbtFileMerge = NULL, hapDir = NULL)
}
\arguments{
  \item{projectName}{
Name given to GBS pipeline. A folder of this name containing folders called "tbt" and "tagcounts" will be created in your working directory if it doesn't already exist. All output from the GBS pipeline will be written in these folders.
}
  \item{keyFile}{
The path to your key file.
}
  \item{seqDir}{
The path to your sequence data.
}

  \item{tasselJar}{
Path to sTASSEL.jar from Tassel 4.0 standalone or an equivalent .jar file compiled from Tassel 4.0 code.
}

  \item{heapSize}{
Maximum heap size for the JVM in GB. Should be 16 or more.
}

  \item{enzyme}{
The enzymes used for restriction digest. Only "PstI-MspI" has been tested.
}
  \item{maxMissingData}{
The maximum allowable percentage of lines with missing values in a returned SNP.
}
  \item{minorAlleleFreq}{
The minimum allowable minor allele frequency for a returned SNP.
}
  \item{maxHet}{
The maximum allowable percentage of heterozygotes for a returned SNP.
}
  \item{isDHpopulation}{
Are the lines from a double haploid population. If yes, heterozygous SNPs are rejected.
}
  \item{callHets}{
Should heterozygous SNPs be returned. If false, heterozygous values for SNPs are returned as missing.
}
  \item{pVal}{
Significance value used for Fisher's test to determine if a SNP should be returned.
}
  \item{runPipeline}{
Should the pipeline be run. If false, a Java object of GBSPipeline class is returned. See details.
}
  \item{isQseq}{
Is sequence data in Qseq format. If false, Fastq format is assumed. Not used if runPipeline is false.
}
  \item{tagCountsDir}{
Directory where tag count files are stored. If NULL, the function file.path(getwd(), projectName, "tagcounts") will be used to set this value.
}
  \item{masterTagsFile}{
Name of the master tags file. If NULL, the function file.path(getwd(), projectName, paste("MasterTags_", projectName, ".cnt", sep="")) will be used to set this value.
}
  \item{tbtDir}{
The name of the tbt directory. If NULL, the function file.path(getwd(), projectName, "tbt") will be used to set this value.
}
  \item{tbtFile}{
Name of the tbt file. If NULL, the function file.path(getwd(), projectName, paste("tbt_", projectName, "_", getDate(), ".bin", sep="")) will be used to set this value.
}
  \item{tbtFileMerge}{
Name of the tbtMerge file. If NULL, the function file.path(getwd(), projectName, paste("tbtMerge_", projectName, "_", getDate(), ".bin", sep="")) will be used to set this value.
}
  \item{hapDir}{
The directory where ".hap" files will be written. If NULL, the function file.path(getwd(), projectName) will be used to set this value.
}
}
\details{
In your current working directory a folder is created with the same name as projectName if it doesn't already exist. Within that folder, folders called "tbt" and "tagcounts" are created if they don't already exist. The arguments are then used to create a Java object of class GBSPipeline using the 'rJava' function \code{\link[rJava]{.jnew}}. If runPipeline=F, this object is returned and methods can manually be performed on this object using the 'rJava' function \code{\link[rJava]{.jcall}}. If runPipeline=T, a sequence of Java methods are performed on your GBSPipeline class to carry out the GBS pipeline. In order, these methods are: "printArgs", "runFastqToTagCountPlugin" (if isQseq=T, "runQseqToTagCountPlugin" is performed instead), "runMergeMultipleTagCountPlugin", "runFastqToTBTPlugin" ("runQseqToTBTPlugin" if isQseq=T), "runMergeTagsByTaxaFilesPlugin", "mergeTaxaInTBT", and "runTagsToSNPsNoAnchor".

The SNP data is outputted in three files in yourProject directory with the names: [yourProject]1_[date].hap, [yourProject]2_[date].hap, and [yourProject]3_[date].hap. The numbers 1, 2 and 3 in the ".hap" files correspond to the number of allowable mismatches during sequence alignment. The three files share many of the same SNPs, but each can contain SNPs not present in either of the other files.
}
\section{Warning}{
Heap size can only be set once per R session, so you must start a new R session if the JVM was launched with an insufficient heap size.

Tassel 4.0 requires Java 1.7 and won't run with other versions of Java, so this version of Java is required to run GBSPipeline.
}
\references{
\cite{Bradbury, P.J., Z. Zhang, D.E. Kroon, T.M. Casstevens, Y. Ramdoss, and E.S. Buckler. 2007. TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635.}

\cite{Poland, J.A., P.J. Brown, M.E. Sorrells, and J.-L. Jannink. 2012. Development of High-Density Genetic Maps for Barley and Wheat Using a Novel Two-Enzyme Genotyping-by-Sequencing Approach. PLoS ONE 7(2): e32253.}
}
\examples{
\dontrun{
#Data to run this example is not included
library(GSwGBS)
setwd("/where/you/want/your/output/directory")
GBSPipeline(projectName="nameOfOutputDirectory",
            keyFile="/path/to/keyfile.txt",
            seqDir="/path/to/sequence/data",
            tasselJar="/path/to/your/tassel4-standalone/sTASSEL.jar", 
            heapSize=48) #heapSize in number of GB Tassel can use
}
}