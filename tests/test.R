library(miRACLe)
#####################################################################################################

   ### Part II: DATA INPUT ###

#The input data contains two parts.
#The first part is the sequence-based information obtained from Targetscan, which is CWCS
seqScore = as.matrix(read.table("./data/conserved_CWCS.txt", head = TRUE, sep = "\t"))

#The other part is expression data
#mirExpr and tarExpr are the mirna and mrna expression data, whose format should be strictly satisfied
#mirtarList provides the relationship of sample names in mirna expression and sample names in mrna expression
mirtarList = as.matrix(read.table("./data/test_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))	#测试数据
mirExpr = as.matrix(read.table("./data/test_miRNA_expression.txt", head = FALSE, sep = "\t"))	#测试数据
tarExpr = as.matrix(read.table("./data/test_mRNA_expression.txt", head = FALSE, sep = "\t"))	#测试数据

mirSelect = "hsa-let-7a-5p"
tarSelect = "PTEN"
samSelect = c("TCGA-FA-A4BB-01A-11R-A31S-13", "TCGA-FA-A4XK-01A-11R-A31S-13", "TCGA-FA-A6HN-01A-11R-A31S-13")

multi_mirSelect = as.matrix(read.table("./data/miRNA_list.txt", head = FALSE, sep = "\t"))
multi_tarSelect = as.matrix(read.table("./data/target_list.txt", head = FALSE, sep = "\t"))

multi_mirSelect = c("hsa-let-7a-5p", "hsa-miR-29a-3p", "hsa-miR-200a-3p")
multi_tarSelect = c("NRAS", "CD274", "CCND2")
#####################################################################################################

   ### Part III: MAIN PROGRAM ###

#主程序1，计算所有miRNA-target的预测分数及排名，基本输入为seqScore是序列匹配评分，mirtarList是miRNA和mRNA表达谱文件中样本的对应关系，mirExpr是miRNA表达谱，tarExpr是mRNA表达谱。
#另有两个参数，exprFilter是表达谱筛选，即在x%的样本中没有表达的miRNA和mRNA会被去掉，默认=1，即不筛选。若选择=0.8，即80%的样本中不表达的miRNA和mRNA将会被去掉。
#ID_pos是样本筛选，默认不筛选，如果用户只想看某一些样本的结果，那么可以通过输入样本id来实现，例如第549行输入了三个样本的id，那么程序就会分析这个三个样本相关的结果，包含individual-level和populationi-level的结果。

final_output = miracle(seqScore, mirtarList, mirExpr, tarExpr)	#默认
final_output = miracle(seqScore, mirtarList, mirExpr, tarExpr, exprFilter=1, ID_pos = samSelect)

write.table(final_output[[2]],file = paste("Individual-level predictions.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)	#输出Individual-level结果
write.table(final_output[[1]],file = paste("Population-level predictions.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)	#输出Population-level结果

#主程序2，计算某个（单个）miRNA对应的所有的靶mRNA分子。基本输入为mirSelect是感兴趣的miRNA的id（如第547行），seqScore是序列匹配评分，mirtarList是miRNA和mRNA表达谱文件中样本的对应关系，mirExpr是miRNA表达谱，tarExpr是mRNA表达谱。
#另有两个参数，exprFilter是表达谱筛选，即在x%的样本中没有表达的miRNA和mRNA会被去掉，默认=1，即不筛选。若选择=0.8，即80%的样本中不表达的miRNA和mRNA将会被去掉。
#ID_pos是样本筛选，默认不筛选，如果用户只想看某一些样本的结果，那么可以通过输入样本id来实现，例如第549行输入了三个样本的id，那么程序就会分析这个三个样本相关的结果，包含individual-level和populationi-level的结果。

mirna_prediction = miracle_mirna(mirSelect, seqScore, mirtarList, mirExpr, tarExpr)	#默认
mirna_prediction = miracle_mirna(mirSelect, seqScore, mirtarList, mirExpr, tarExpr, exprFilter=0.9, ID_pos = samSelect)

write.table(mirna_prediction[[1]],file = paste("Individual-level prediction result of ", mirSelect, ".txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)	#输出Individual-level结果
write.table(mirna_prediction[[2]],file = paste("Population-level prediction result of ", mirSelect, ".txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)	#输出Population-level结果

#主程序3，计算某个（单个）mRNA对应的所有的miRNA分子。基本输入为tarSelect是感兴趣的mRNA（如第548行），seqScore是序列匹配评分，mirtarList是miRNA和mRNA表达谱文件中样本的对应关系，mirExpr是miRNA表达谱，tarExpr是mRNA表达谱。
#另有两个参数，exprFilter是表达谱筛选，即在x%的样本中没有表达的miRNA和mRNA会被去掉，默认=1，即不筛选。若选择=0.8，即80%的样本中不表达的miRNA和mRNA将会被去掉。
#ID_pos是样本筛选，默认不筛选，如果用户只想看某一些样本的结果，那么可以通过输入样本id来实现，例如第549行输入了三个样本的id，那么程序就会分析这个三个样本相关的结果，包含individual-level和populationi-level的结果。

mrna_prediction = miracle_target(tarSelect, seqScore, mirtarList, mirExpr, tarExpr)	#默认
mrna_prediction = miracle_target(tarSelect, seqScore, mirtarList, mirExpr, tarExpr, exprFilter=0.8, ID_pos = samSelect)

write.table(mrna_prediction[[1]],file = paste("Individual-level prediction result of ", tarSelect, ".txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)	#输出Individual-level结果
write.table(mrna_prediction[[2]],file = paste("Population-level prediction result of ", tarSelect, ".txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)	#输出Population-level结果

#主程序4，计算某些（多个）miRNA对应的所有的靶mRNA分子。基本输入为multi_mirSelect是感兴趣的miRNA的id list（如第554行），seqScore是序列匹配评分，mirtarList是miRNA和mRNA表达谱文件中样本的对应关系，mirExpr是miRNA表达谱，tarExpr是mRNA表达谱。
#另有两个参数，exprFilter是表达谱筛选，即在x%的样本中没有表达的miRNA和mRNA会被去掉，默认=1，即不筛选。若选择=0.8，即80%的样本中不表达的miRNA和mRNA将会被去掉。
#ID_pos是样本筛选，默认不筛选，如果用户只想看某一些样本的结果，那么可以通过输入样本id来实现，例如第549行输入了三个样本的id，那么程序就会分析这个三个样本相关的结果，包含individual-level和populationi-level的结果。

multi_mirna_prediction = miracle_multi_mirna(multi_mirSelect, seqScore, mirtarList, mirExpr, tarExpr)	#默认
multi_mirna_prediction = miracle_multi_mirna(multi_mirSelect, seqScore, mirtarList, mirExpr, tarExpr, exprFilter=0.8, ID_pos = samSelect)

write.table(multi_mirna_prediction[[1]],file = paste("Individual-level prediction result of selected miRNAs.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(multi_mirna_prediction[[2]],file = paste("Population-level prediction result of selected miRNAs.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)

#主程序5，计算某些（多个）mRNA对应的所有的miRNA分子。基本输入为multi_tarSelect是感兴趣的mRNA的symbol list（如第555行），seqScore是序列匹配评分，mirtarList是miRNA和mRNA表达谱文件中样本的对应关系，mirExpr是miRNA表达谱，tarExpr是mRNA表达谱。
#另有两个参数，exprFilter是表达谱筛选，即在x%的样本中没有表达的miRNA和mRNA会被去掉，默认=1，即不筛选。若选择=0.8，即80%的样本中不表达的miRNA和mRNA将会被去掉。
#ID_pos是样本筛选，默认不筛选，如果用户只想看某一些样本的结果，那么可以通过输入样本id来实现，例如第549行输入了三个样本的id，那么程序就会分析这个三个样本相关的结果，包含individual-level和populationi-level的结果。

multi_mrna_prediction = miracle_multi_target(multi_tarSelect, seqScore, mirtarList, mirExpr, tarExpr)	#默认
multi_mrna_prediction = miracle_multi_target(multi_tarSelect, seqScore, mirtarList, mirExpr, tarExpr, exprFilter=0.9, ID_pos = samSelect)

write.table(multi_mrna_prediction[[1]],file = paste("Individual-level prediction result of selected targets.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(multi_mrna_prediction[[2]],file = paste("Population-level prediction result of selected targets.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)

