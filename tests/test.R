#####################################################################################################

### Part II: DATA INPUT ###

#The input data contains two parts.
#The first part is the sequence-based information obtained from Targetscan, which is CWCS
seqScore = as.matrix(read.table("CWCS.txt", head = TRUE, sep = "\t"))

#The other part is expression data
#mirExpr and taExpr are the mirna and mrna expression data, whose format should be strictly satisfied
#mirtaList provides the relationship of sample names in mirna expression and sample names in mrna expression
mirtaList = as.matrix(read.table("DLBC_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))
mirExpr = as.matrix(read.table("DLBC_miRNA_expression.txt", head = FALSE, sep = "\t"))
taExpr = as.matrix(read.table("DLBC_mRNA_expression.txt", head = FALSE, sep = "\t"))
#topList is the number of top rankers that users are interested in.
topList = 5000

#####################################################################################################

### Part III: MAIN PROGRAM ###

#Our algorithm provides both population-level result and sample-level result
final_output = miracle_score(seqScore, mirtaList, mirExpr, taExpr, topList)

write.table(final_output[[1]],file = paste("Population-level top", topList, "predictions.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_output[[2]],file = "Population-level top 1 percent predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_output[[3]],file = paste("Individual-level top", topList, "predictions.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_output[[4]],file = "Individual-level top 1 percent predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
