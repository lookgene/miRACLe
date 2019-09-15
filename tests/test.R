library(miRACLe)
### Part II: DATA INPUT ###

# The input data contains two parts.
# The first part is the sequence-based interaction scores obtained from sequence-based predictors, default
# is cumulative weighted context++ score for conserved target sites of conserved miRNA families obtained
# from TargetScan v7.2.

seqScore = as.matrix(read.table("../TargetScan7_CWCS_cons.txt", head = TRUE, sep = "\t"))

# The other part is expression data
# mirExpr and tarExpr refer to the mirna and mrna expression data, respectively
# sampleMatch provides the corresponding relationship of sample names in mirna expression and sample names in mrna expression
sampleMatch = as.matrix(read.table("Test_sampleMatch.txt", head = TRUE, sep = "\t"))	#test data
mirExpr = as.matrix(read.table("Test_miRNA_expression.txt", head = FALSE, sep = "\t"))	#test data
tarExpr = as.matrix(read.table("Test_mRNA_expression.txt", head = FALSE, sep = "\t"))	#test data

sampleSelect = c("TCGA-FA-A4BB-01A-11R-A31S-13", "TCGA-FA-A4XK-01A-11R-A31S-13", "TCGA-FA-A6HN-01A-11R-A31S-13")

mirExpr_ind = as.matrix(read.table("Test_ind_miRNA.txt", head = FALSE, sep = "\t"))	#test data
tarExpr_ind = as.matrix(read.table("Test_ind_mRNA.txt", head = FALSE, sep = "\t"))	#test data

#####################################################################################################

### Part III: MAIN PROGRAM ###

# function 1: calculate the predicted score and rank of all miRNA-target pairs, the essential inputs are seqScore
# (sequence-based interaction scores), sampleMatch(corresponding relationships between samples from miRNA expression
# data and mRNA expression data), mirExpr(the expression profile of miRNA), tarExpr(the expression profile of mRNA),
# another three optional inputs are: exprFilter(filter of expression profile, default is 1, miRNAs/mRNAs that are not
# expressed in more than a given percentage of samples will be removed), samSelect(selection of samples, users can
# select a part of all samples to analyze, default is no selection applied), OutputSelect(logical, select “TRUE”
# to return the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result.
# Default is TRUE.)

final_output = miracle(seqScore, sampleMatch, mirExpr, tarExpr)	#default
final_output = miracle(seqScore, sampleMatch, mirExpr, tarExpr, exprFilter = 1, samSelect = sampleSelect, OutputSelect = TRUE)

write.table(final_output$Ind, file = "Individual-level predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)	#Individual-level result
write.table(final_output$Pop, file = "Population-level predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)	#Population-level result

# function 2: calculate the predicted score and rank of all miRNA-target pairs for individual sample, the essential
# inputs are seqScore(sequence-based interaction scores), mirExpr(the expression profile of miRNA), tarExpr
# (the expression profile of mRNA), another one optional input is OutputSelect(logical, select “TRUE” to return
# the top 10 percent-ranked predictions by scores, and “FALSE” to return the whole prediction result.
# Default is TRUE.)

final_output_ind = miracle_ind(seqScore, mirExpr_ind, tarExpr_ind)	#default
final_output_ind = miracle_ind(seqScore, mirExpr_ind, tarExpr_ind, OutputSelect = TRUE)

write.table(final_output_ind, file = "prediction result of one sample.txt", quote = FALSE, sep = "\t", row.names = FALSE)
