# Tutorial Code

#Procuring Data
# NCBI enables batch downloads, in which case all sequences will be in a single fasta. In the case of GISAID, likely many files will need to be downloaded and then merged via unix commands
# create a directory with the fasta files and move all the fasta into there. Navigate to that directory in unix and write this command to compile them all into one file.

# now these are ready for alignment

#### ALIGNMENT STEP ###
# BiocManager::install("DECIPHER")
library(DECIPHER)

fas <- "merged.fasta"
fas <- "NAclades.fasta"
fas <- "supertree_accessions.fasta"
fas <- "full_supertree_accessions.fasta"
seqs <- readDNAStringSet(fas)
aligned <- AlignSeqs(seqs)

### Write to Fasta####

# fasta need to be written in a certain format for it to be interpreted by the PCA analysis well.

# creating a custom function to do so. 
alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  n <- length(names(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', names(alignment)[i], "\t"))
    the.sequence <- as.character(aligned[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}
# function takes in the alignment object and the desired output file name 

alignment2Fasta(seqs, 'full_supertree_aligned.fasta')

### Read in fasta as tab-delimited file###
sites <- read.table(file="full_supertree_aligned.fasta", header=F, sep="\t")
sites <- read.table(file="PCA_FINAL_NA.fasta", header=F, sep="\t")
dim(sites)
# we have a matrix where the header is in the first column, sequence is in the second column

### finding the size of data
n.sample <- dim(sites)[1] # number of samples
seq.len <- nchar(sites[2,2]) # length of the sequences

# translating the vector into a boolean array
boolean.matrix <- array(0, dim=c(n.sample, 5*seq.len))

colnames(boolean.matrix) <- c(paste("A_", 1:seq.len, sep=""),paste("T_", 1:seq.len, sep=""),paste("G_", 1:seq.len, sep=""),paste("C_", 1:seq.len, sep=""),paste("N_", 1:seq.len, sep=""))

# setting the row names of the matrix to hold the respective header information
rownames(boolean.matrix) <- sites[ ,1]


# populating the matrix with one-hot encoding 
for (s in 1:n.sample){
  se <- sites[s, 2]
  se <- tolower(se)
  
  for (le in  1:seq.len){
    base <- substr(se, le,le)
    
    if(base =="a") {
      boolean.matrix[s, le] <-1
    } else {
      
      if(base =="t") {
        boolean.matrix[s, le+seq.len] <-1
      } else {
        
        if(base =="g") {
          boolean.matrix[s, le+seq.len*2] <-1
        } else {
          
          if(base =="c") {
            boolean.matrix[s, le+seq.len*3] <-1
          } else {
            
            boolean.matrix[s, le+seq.len*4] <-1
          }}}}
  }}

# checking that all nucleotides are accounted for.
apply(boolean.matrix, 1, sum)

###PRINCIPAL COMPONENT ANALYSIS####
##################################
## centering the data for PCA.
center<- apply(boolean.matrix, 2, mean) # extracting the mean of all columns
# scaling the vectors (centering) and creating a distance matrix
diffs<-sweep(boolean.matrix, 2, center)

# compensating for the doubled counts in Euclidean distance metrics
diffs <- diffs/(2^0.5)   

# checking distribution of the distances
distances<- (apply(diffs^2, 1, sum))^0.5
qqnorm(distances)

library(ggplot2)

# histogram of distances
ggplot(data=as.data.frame(distances), aes(x=distances))+
  geom_histogram(color="darkblue", fill="lightblue", bins = 20)

#### PCA core
# singular value decomposition of a matrix (similar to identifying the linear combinations in a matrix)
res_svd <- svd(diffs)  #
# d: vector of singular values, u: matrix with columns of left singular values, v: matrix with columns of right singular values.
str(res_svd)
Left <- res_svd$u		# the left singular vector
Right <- res_svd$v		# the right singular vector
sqL <- diag(res_svd$d)		# diagonal matrix of the singular values (identifying the variance from the principal components)

### calculation of pc's 
sPC_nuc  	<-	 Right %*% sqL / (n.sample^0.5)
sPC_sample	 <-	 Left %*% sqL/ (seq.len^0.5)

rownames(sPC_nuc)<- colnames(boolean.matrix) 
rownames(sPC_sample)<- rownames(boolean.matrix) 

#### output to text files for documentation or manual inspection. 
write.table(sPC_sample, file="sPC_sample.txt", sep="\t")
write.table(sPC_nuc, file="sPC_nuc.txt", sep="\t")

write.table(sPC_sample, file="full_supertree_sPC_sample.txt", sep="\t")
write.table(sPC_nuc, file="full_supertree_sPC_nuc.txt", sep="\t")

### VISUALIZATIONS###
# scree plot (contribution of each principal component towards explaining the overall variation in the data)

index <- 1:20
contribution.scores <- round((res_svd$d/sum(res_svd$d)*100)[1:20],2)
contributions <- data.frame(index,contribution.scores)

ggplot(data=contributions, aes(x=`index`, y=`contribution.scores`, group=1)) +
  geom_line(linetype = "twodash")+
  geom_text(aes(label=`contribution.scores`),hjust=-0.54, vjust=-0.6, size=3)+
  geom_point() + labs(y="PC Contribution (%)", x = "PC") + ggtitle("Cumulative Proportion of Variance Described by Leading 20 PCs") 

#### PCA plot of the samples
ggplot(as.data.frame(sPC_sample), aes(x=V1, y=V2)) +
  geom_point(size=2) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples") 

sample.graham <- read.table("full_supertree_sPC_sample.txt", sep="\t")
rownames((sample.graham))
names(sample.graham)
sample.graham$accessions <- rownames((sample.graham))
sample.graham$accessions <- str_extract(sample.graham$accessions,"[A-Z]+[0-9]+")

ggplot(as.data.frame(sample.graham), aes(x=V1, y=V2, label=`accessions`)) +
  geom_text(size=3) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples")

type <- c("SARS","SARS","SARS","SARS","BAT","BAT", "BAT", "SARS", "SARS", "MERS", "BAT", "BAT", "BAT", "MERS", "BAT", "BAT", "BAT", "BAT", "SARS-COV-2", "REF", rep("SARS-COV-2", 100))

sample.graham$type <- as.factor(type)
class(sample.graham$type)

ggplot(data=sample.graham, aes(x=V1, y=V2, color=`type`)) +
  geom_point(size=2) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples") 


colors <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4")
colors <- colors[as.numeric(sample.graham$type)]

s3d <- scatterplot3d(sample.graham[1:3], pch=19, main="3D Principal Component Plot of SARS-CoV2 Clades", color=colors)
legend("right", legend = levels(sample.graham$type),
       col =  c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4"), pch = 16)

# PCA plot of the positions
position.index <- 1:seq.len
A.<-sPC_nuc[1:seq.len,1]
T.<-sPC_nuc[1:seq.len+seq.len,1]
G. <- sPC_nuc[1:seq.len+seq.len+seq.len,1]
C. <- sPC_nuc[1:seq.len+seq.len+seq.len+seq.len,1]

sites.plotting <- data.frame(position.index, A., G., C., T.)
library(reshape2)
molten.data<-melt(sites.plotting,
                  measure.vars=c("A.", "G.", "C.", "T."),
                  variable_name = "variable")
molten.data$label <- paste(molten.data$variable, "_", molten.data$position.index)

ggplot(molten.data, aes(x=`position.index`, y=`value`)) + geom_point(aes(color = `variable`, shape = `variable`), size = 2) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) + theme(legend.position="bottom")+
  labs(y="sPC1", x = "Sites") + ggtitle("Site Specific Variation") 

# or use label names to see the exact position numbers.
ggplot(molten.data, aes(x=`position.index`, y=`value`)) + geom_text(aes(color = `variable`, label=`label`), size = 1) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) + theme(legend.position="bottom")+
  labs(y="sPC1", x = "Sites") + ggtitle("Site Specific Variation") 

cat *.fasta > combined.fasta
# extracting isolate names
test <- rownames(sPC_sample)
print(test)
testing <- gsub("^.*-2/","", test)
print(testing)

testing <- gsub("2020,.*","", testing)

str_extract(label.testing,"[A-Z]+[0-9]+")

sample.pcs <- data.frame(sPC_sample, labelling)
ggplot(as.data.frame(sample.pcs), aes(x=X1, y=X2, label=`labelling`)) +
  geom_text(size=3) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples") 
library(stringr)
accessions <- str_extract(labels,"[A-Z]+[0-9]+")

sample.pcs <- data.frame(sPC_sample, accessions)

ggplot(as.data.frame(sample.pcs), aes(x=X1, y=X2, label=`accessions`)) +
  geom_text(size=3) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples") 

test <- ">hCoV-19/Germany/HH-HPI-p4856/2021|EPI_ISL_1464582|2021-03-04"
str_extract(test,"EPI_ISL_[0-9]+")
sample.pcs <- data.frame(sPC_sample, accessions)
accessions <- rownames(sPC_sample)

accessions <- str_extract(accessions,"[A-Z]+[0-9]+")

ggplot(as.data.frame(sample.pcs), aes(x=X1, y=X2, label=`accessions`)) +
  geom_text(size=3) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples")
ggplot(as.data.frame(sample.pcs), aes(x=X2, y=X3, label=`accessions`)) +
  geom_text(size=3) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples")
ggplot(as.data.frame(sample.pcs), aes(x=X2, y=X3, label=`accessions`)) +
  geom_text(size=3) + labs(y="sPC2", x = "sPC1") + ggtitle("PCs of Samples")


sample.pcs$clade <- c("L", "GRY", "GRY", "GRY", "GR", "GR", "GR", "G", "G", "G", "S", "O", "L", "S", "S", "O", "GH", "GH", "O", "V", "GV", "V", "GV", "GV", "S", "S", "S", "GH", "V", "L")

library(scatterplot3d)
scatter <- read.table("sample.pcs.txt", sep = "\t")
scatterplot3d(scatter[1:3], pch=19, color = as.numeric(scatter$clade), main="3D Principal Component Plot of SARS-CoV2 Clades")
legend("right", legend = levels(scatter$clade), col =  c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352", "#CC79A7"), pch = 16)



