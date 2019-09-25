#' miRACLe predicts miRNA-mRNA interactions for individual samples as well as
#' sample population
#'
#' Infer miRNA-mRNA interactions using paired miRNA-mRNA expression profile.
#' miRACLe is based on a random contact model. It integrates genome-wide
#' expression profiles and sequence-based interaction scores of miRNAs and mRNAs
#' to calculate the relative probability of effective contacts. We consider
#' that the more likely an effective contact takes place between two molecules,
#' the more likely that a regulatory relationship between them exists.
#' This function supports prediction of miRNA-mRNA interactions at both
#' individual and population levels.
#'
#' @param seqScore  A Z x 3 data.frame that contains sequence-based interaction
#' scores for putative miRNA-mRNA pairs. The first column must contain gene
#' symbols for mRNAs, the second column must contain miRNA ids for miRNAs, the
#' third column must contain sequence-based score for the miRNA-mRNA pair
#' @param sampleMatch A character K x 2 matrix of paired sample identifiers used
#'  in mirExpr and tarExpr. The first column must contain the sample identifiers
#'   used in mirExpr and the second column must contain the sample identifiers
#'   used in tarExpr. Sample identifiers in the same row should be paired.
#' @param mirExpr A numeric M x R matrix as the observed expression of M miRNAs
#' in R samples. NOTE: colnames is required for mirExpr to map each sample to
#' the sample identifiers in the first column of sampleMatch. rownames is
#' required for mirExpr to map each miRNA to the miRNA id in the first column
#' of sampleMatch.
#' @param tarExpr A numeric N x T matrix as the observed expression of N mRNAs
#' in T samples. NOTE: colnames is required for tarExpr to map each sample to
#' the sample identifiers in the second column of sampleMatch. rownames is
#' required for tarExpr to map each mRNA to the gene symbol in the first column
#'  of sampleMatch.
#' @param samSelect Optional, a vector of samples that users select to get the
#' integrated result over those samples. The elements in this vector should be
#' matched to the sample identifiers in mirExpr. Default is NULL, which means
#' that all samples in sampleMatch will be analyzed.
#' @param exprFilter a filter to remove the lowly expressed miRNAs and mRNAs.
#' If provided, exprFiter=k, then both miRNAs and mRNAs that are expressed in
#' less than (1-k)% of the samples in the provided datasets will be filtered
#' out. Default is 1.
#' @param OutputSelect Logical, select “TRUE” to return the top 10 percent
#' ranked predictions by scores, and “FALSE” to return the whole prediction
#' result. Default is TRUE.
#'
#' @section Detail:
#' miRACLe predicts miRNA-mRNA interactions for individual samples as well as
#' sample population. For individual-level prediction, miRACLe inputs a single
#' paired miRNA-mRNA expression profile and sequence-based scores for putative
#' miRNA-mRNA pairs and output a list of miRNA-mRNA interactions ranked by
#' miRACLe scores. When a sample population is of interest, the individualized
#' miRACLe scores will be averaged over all samples in the population to
#' determine a population-level miRACLe score, based on which we output the
#' ranking list of miRNA-mRNA interactions.
#'
#' In order to run the current version of miRACLe, the users should provide two
#' data files that describe the expression levels of each miRNA and mRNA. And
#' one additional file that defines the correspondence of samples between the
#' miRNA and mRNA data files. All files are tab-delimited ASCII text files. As
#' for expression data, both microarray profiling and RNA sequencing data are
#' accepted. To achieve optimal prediction on the sequencing data, log2
#' transformed normalized read counts (e.g. RSEM) are recommended as the input
#' for the program. In order for this function to execute correctly, this
#' expression matrix should be transformed into a non-negative numeric matrix.
#'
#' miRACLe provides sequence-based interaction scores for putative miRNA-mRNA
#' pairs. These scores are originally obtained from TargetScan v7.2
#' (TargetScan_CWCS_cons and TargetScan_CWCS), DIANA-microT-CDS (DIANA_microT_CDS),
#' MirTarget v4 (MirTarget4), miRanda-mirSVR (miRanda_mirSVR) and compiled by the
#' developers to fit the model. Default is TargetScan_CWCS_cons. The other scores
#' can be downloaded \href{https://figshare.com/s/0b7c68cd5152da27a191}{here}.
#' User can also provide their own seqScore, as long as the format meets the
#' requirements.
#'
#' @return An object list containing the following items:
#'
#' Ind: inferred individual-level prediction result, a 4*N-column
#'
#' matrixPop: inferred population-level prediction result, a 4-column matrix
#'
#' @export
#' @section Note:
#' Sample identifiers in sampleMatch should be a subset of those in mirExpr and tarExpr.
#' @examples
#' data(seqScore)   # load the default sequence-based interaction score
#' data(Test_data)  # load test datasets
#' mirExpr <- Test_DLBC_miRNA # miRNA expression
#' tarExpr <- Test_DLBC_mRNA  # mRNA expression
#' sampleMatch <- Test_DLBC_sampleMatch # sample matching file
#' sampleSelect = c("TCGA-FA-A4BB-01A-11R-A31S-13", "TCGA-FA-A4XK-01A-11R-A31S-13", "TCGA-FA-A6HN-01A-11R-A31S-13") # samples selected from the test dataset to analyze
#' final_output <- miracle(seqScore, sampleMatch, mirExpr, tarExpr, samSelect = sampleSelect, exprFilter = 0.8, OutputSelect = FALSE)
#' final_output$Ind   # Individual-level result
#' final_output$Pop   # Population-level result
miracle = function(seqScore, sampleMatch, mirExpr, tarExpr,samSelect = NULL, exprFilter = 1, OutputSelect = TRUE){
  input_list = matrix_transfer(seqScore)
  x = name_func(sampleMatch, mirExpr, tarExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(tarExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], exprFilter)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  if(length(samSelect) == 0){
    z9 = miracle_detail(z3[[1]], z3[[2]], z7[[3]], z8, sampleMatch[, 1], z3[[3]], z3[[4]], OutputSelect = TRUE, select_ID = sampleMatch[, 1])
  }
  else{
    z9 = miracle_detail(z3[[1]], z3[[2]], z7[[3]], z8, sampleMatch[, 1], z3[[3]], z3[[4]], OutputSelect = TRUE, select_ID = samSelect)
  }
  names(z9) = c("Ind", "Pop")
  return(z9)
}

#' miRACLe predicts miRNA-mRNA interactions for individual samples
#'
#' Infer miRNA-mRNA interactions using a single paired miRNA-mRNA expression profile.
#' miRACLe is based on a random contact model. It integrates genome-wide expression
#' profile and sequence-based interaction scores of miRNAs and mRNAs to calculate
#' the relative probability of effective contacts. We consider that the more likely
#' an effective contact takes place between two molecules, the more likely that a
#' regulatory relationship between them exists. This function is specifically
#' designed to predict individualized miRNA-mRNA interactions using a single
#' paired miRNA-mRNA expression profile.
#'
#' @param seqScore A Z x 3 data.frame that contains sequence-based interaction
#' scores for putative miRNA-mRNA pairs. The first column must contain gene
#' symbols for mRNAs, the second column must contain miRNA ids for miRNAs,
#' the third column must contain sequence-based score for the miRNA-mRNA pair.
#' @param mirExpr A numeric vector as the observed expression of M miRNAs. NOTE:
#'  rownames is required for mirExpr_ind to map each miRNA to the miRNA id in
#'  the second column of seqScore.
#' @param tarExpr A numeric vector as the observed expression of N mRNAs. NOTE:
#'  rownames is required for tarExpr_ind to map each mRNA to the gene symbol in
#'  the first column of seqScore.
#' @param OutputSelect Logical, select “TRUE” to return the top 10 percent
#' ranked predictions by scores, and “FALSE” to return the whole prediction
#' result. Default is TRUE
#'
#' @return An object list containing the following items:
#'
#' Ind: inferred prediction result, a 4-column matrix
#' @section Detail:
#' miRACLe predicts miRNA-mRNA interactions for individual samples as well as
#' sample population. For individual-level prediction, miRACLe inputs a single
#' paired miRNA-mRNA expression profile and sequence-based scores for putative
#' miRNA-mRNA pairs and output a list of miRNA-mRNA interactions ranked by
#' miRACLe scores. When a sample population is of interest, the individualized
#' miRACLe scores will be averaged over all samples in the population to
#' determine a population-level miRACLe score, based on which we output the
#' ranking list of miRNA-mRNA interactions.
#'
#' In order to run the current version of miRACLe, the users should provide two
#' data files that describe the expression levels of each miRNA and mRNA. And
#' one additional file that defines the correspondence of samples between the
#' miRNA and mRNA data files. All files are tab-delimited ASCII text files. As
#' for expression data, both microarray profiling and RNA sequencing data are
#' accepted. To achieve optimal prediction on the sequencing data, log2
#' transformed normalized read counts (e.g. RSEM) are recommended as the input
#' for the program. In order for this function to execute correctly, this
#' expression matrix should be transformed into a non-negative numeric matrix.
#'
#' miRACLe provides sequence-based interaction scores for putative miRNA-mRNA
#' pairs. These scores are originally obtained from TargetScan v7.2
#' (TargetScan_CWCS_cons and TargetScan_CWCS), DIANA-microT-CDS (DIANA_microT_CDS),
#' MirTarget v4 (MirTarget4), miRanda-mirSVR (miRanda_mirSVR) and compiled by the
#' developers to fit the model. Default is TargetScan_CWCS_cons. The other scores
#' can be downloaded \href{https://figshare.com/s/0b7c68cd5152da27a191}{here}.
#' User can also provide their own seqScore, as long as the format meets the
#' requirements.
#' @export
#'
#' @examples
#' data(seqScore)   # load the default sequence-based interaction score
#' data(Test_data)  # load test datasets
#' mirExpr_ind <- Test_Hela_miRNA # miRNA expression
#' tarExpr_ind <- Test_Hela_mRNA  # mRNA expression
#' final_output_ind <- miracle_ind(seqScore, mirExpr_ind, tarExpr_ind, OutputSelect = TRUE)
miracle_ind = function(seqScore, mirExpr, tarExpr, OutputSelect = TRUE){
  input_list = matrix_transfer(seqScore)
  z1 = mirna_matrix_ind(mirExpr, input_list[[1]])
  z2 = mrna_matrix_ind(tarExpr, input_list[[2]])
  z3 = mirna_mrna_data_ind(z1[[1]], z2[[1]], z1[[2]], z2[[2]])
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  Ind = miracle_detail_ind(z3[[1]], z3[[2]], z7[[3]], z8, z3[[3]], z3[[4]], OutputSelect)
  return(Ind)
}


matrix_transfer = function(input_matrix){
  mirna_name = unique(input_matrix[, 2])
  mrna_name = unique(input_matrix[, 1])
  mod_matrix = matrix(data = 0, ncol = length(mirna_name), nrow = length(mrna_name))
  pos1 = match(input_matrix[, 2], mirna_name)
  pos2 = match(input_matrix[, 1], mrna_name)
  for(i in 1:nrow(input_matrix)){
    mod_matrix[pos2[i], pos1[i]] = mod_matrix[pos2[i], pos1[i]] + abs(as.numeric(input_matrix[i, 3]))
  }
  return(list(mirna_name, mrna_name, mod_matrix))
}
#The following three functions are used to match samples, mirnas and mrnas
#By using these three functions, we can get the paired mrna-mirna expression data
name_func = function(name_base, mirna_base, mrna_base){
  x = match(name_base[, 1], mirna_base[1, ])
  y = match(name_base[, 2], mrna_base[1, ])
  return(list(x,y))
}
mirna_matrix = function(mirna_base, mirna_name, mirna){
  mirna_name1 = rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i] = unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use = mirna_base[mirna_name1 %in% mirna, mirna_name]
  mirna_name1 = mirna_name1[mirna_name1 %in% mirna]
  return(list(mirna_use, mirna_name1))
}
mirna_matrix_ind = function(mirna_base, mirna){
  mirna_name = rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name[i] = unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use = mirna_base[mirna_name %in% mirna, 2]
  mirna_name = mirna_name[mirna_name %in% mirna]
  return(list(mirna_use, mirna_name))
}
mrna_matrix = function(mrna_base, mrna_name, mrna){
  mrna_exist = rep(0, nrow(mrna_base))
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
  }
  mrna_use = mrna_base[mrna_name1 %in% mrna, mrna_name]
  mrna_name2 = mrna_name1[mrna_name1 %in% mrna]
  return(list(mrna_use, mrna_name2))
}
mrna_matrix_ind = function(mrna_base, mrna){
  mrna_name = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
  }
  mrna_use = mrna_base[mrna_name %in% mrna, 2]
  mrna_name = mrna_name[mrna_name %in% mrna]
  return(list(mrna_use, mrna_name))

}
#mirna_mrna_data will select the RNAs that are expressed in a certain percentage of samples
mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, exprFilter){
  if(is.vector(mirna_use) == TRUE){
    mirna_use = as.matrix(mirna_use)
  }
  if(is.vector(mrna_use) == TRUE){
    mrna_use = as.matrix(mrna_use)
  }
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
  mirna_use[mirna_use < 0] = 0
  mrna_use[mrna_use < 0] = 0
  mirna_sgn = seq(1, nrow(mirna_use), 1)
  mrna_sgn = seq(1, nrow(mrna_use), 1)
  for (i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,] == 0) >= exprFilter * ncol(mirna_use)){
      mirna_sgn[i] = 0
    }
  }
  for (i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,] == 0) >= exprFilter * ncol(mrna_use)){
      mrna_sgn[i] = 0
    }
  }
  mirna_use1 = mirna_use[mirna_sgn, ]
  mrna_use1 = mrna_use[mrna_sgn, ]
  mirna_name1 = mirna_name[mirna_sgn]
  mrna_name1 = mrna_name[mrna_sgn]
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1))
}
mirna_mrna_data_ind = function(mirna_use, mrna_use, mirna_name, mrna_name){
  if(is.vector(mirna_use) == TRUE){
    mirna_use = as.matrix(mirna_use)
  }
  if(is.vector(mrna_use) == TRUE){
    mrna_use = as.matrix(mrna_use)
  }
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
  mirna_use[mirna_use < 0] = 0
  mrna_use[mrna_use < 0] = 0

  return(list(mirna_use, mrna_use, mirna_name, mrna_name))
}
#mirna_mrna_loc is used to get the sequence matrix that fits the expression data
mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, CWCS){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  CWCS = t(CWCS)
  CWCS_use = CWCS[mirna_loc, mrna_loc]
  rownames(CWCS_use) = mirna_name
  colnames(CWCS_use) = mrna_name
  return(list(mirna_loc, mrna_loc, CWCS_use))
}

#the following function sum_miracle and miracle_detail are the core functions of the miracle algorithm
#sum_miracle is used to calculate the normalization coefficient of miracle algorithm
sum_miracle = function(mirna, mrna, CWCS){
  if(is.vector(mirna) == TRUE){
    mirna = as.matrix(mirna)
  }
  if(is.vector(mrna) == TRUE){
    mrna = as.matrix(mrna)
  }

  res = rep(0, ncol(mirna))
  for(i in 1:ncol(mirna)){
    mirna_use1 = as.numeric(mirna[, i])
    mrna_use1 = as.numeric(mrna[, i])
    temp = mirna_use1 %*% t(mrna_use1)
    temp2 = CWCS*temp
    res[i] = sum(temp2)
  }
  return(res)
}
#miracle_basic is used to get the detailed information of the top rankers by miracle score, both population-level and individual-level
#The default is the top 10 percent-ranked pairs, while the users can set the OutputSelect "FALSE" to get the full result
miracle_detail = function(mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, OutputSelect = TRUE, select_ID){
  if(OutputSelect == TRUE){
    thres = 0.1
  }
  else{
    thres = 1
  }
  ID_poss = match(select_ID, ID)
  ID_pos = ID_poss[which(!is.na(ID_poss))]
  if(length(ID_pos) == 0){
    print("The selected samples are not in the predicted list.")
    return("The selected samples are not in the predicted list.")
  }
  else{
    if(is.vector(mirna_use1) == TRUE){
      mirna_use = matrix(as.numeric(mirna_use1), ncol = 1)
    }
    else{
      mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
    }
    if(is.vector(mrna_use1) == TRUE){
      mrna_use = matrix(as.numeric(mrna_use1), ncol = 1)
    }
    else{
      mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
    }
    mirna_use = mirna_use[, ID_pos]
    mrna_use = mrna_use[, ID_pos]
    if(is.vector(mirna_use) == TRUE){
      mirna_use = as.matrix(mirna_use)
    }
    if(is.vector(mrna_use) == TRUE){
      mrna_use = as.matrix(mrna_use)
    }
    numm = floor(length(CWCS[CWCS != 0]) * thres)
    #pro_matrix_agg saves the probability matrix of the integrated result
    pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
    detail_output = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = numm, byrow = TRUE)
    temp_name = c()
    num_max = 0

    for(m in 1:ncol(mirna_use)){
      pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
      pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
      num1 = floor(length(pro_matrix_temp[pro_matrix_temp != 0]) * thres)
      sort_one_rank = order(pro_matrix_temp, decreasing = TRUE)[1 : num1]
      sort_one_prob = sort(pro_matrix_temp, decreasing = TRUE)[1 : num1]
      temp3 = ceiling(sort_one_rank/length(mirna_name))
      temp4 = sort_one_rank - (temp3 - 1) * length(mirna_name)
      for(i in 1 : num1){
        detail_output[i, (4 * m - 3)] = i
        detail_output[i, (4 * m - 2)] = mrna_name[temp3[i]]
        detail_output[i, (4 * m - 1)] = mirna_name[temp4[i]]
        detail_output[i, (4 * m)] = sort_one_prob[i]
      }
      temp_name = c(temp_name, ID[m], "Gene Symbol", "miRNA Name", "miRACLe Score")
      num_max = max(num_max, num1)
    }
    colnames(detail_output) = temp_name
    detail_output1 = detail_output[1:num_max, ]

    pro_matrix_agg = pro_matrix_agg/ncol(mirna_use)
    num2 = floor(length(pro_matrix_agg[pro_matrix_agg != 0]) * thres)
    out_matrix = matrix(data = 0, ncol = 4, nrow = num2, byrow = TRUE)
    temp_rank = order(pro_matrix_agg, decreasing = TRUE)[1 : num2]
    temp_prob = sort(pro_matrix_agg, decreasing = TRUE)[1 : num2]
    mrna_pos = ceiling(temp_rank/length(mirna_name))
    mirna_pos = temp_rank - (mrna_pos - 1) * length(mirna_name)

    for(i in 1:nrow(out_matrix)){
      out_matrix[i, 1] = i
      out_matrix[i, 2] = mrna_name[mrna_pos[i]]
      out_matrix[i, 3] = mirna_name[mirna_pos[i]]
      out_matrix[i, 4] = temp_prob[i]
    }
    colnames(out_matrix) = c("rank", "Gene Symbol", "miRNA Name", "miRACLe Score")

  }
  return(list(detail_output1, out_matrix))
}
miracle_detail_ind = function(mirna_use1, mrna_use1, CWCS, sumup, mirna_name, mrna_name, OutputSelect = TRUE){
  if(OutputSelect == TRUE){
    thres = 0.1
  }
  else{
    thres = 1
  }
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  mirna_use = as.matrix(mirna_use)
  mrna_use = as.matrix(mrna_use)

  pro_matrix = (mirna_use %*% t(mrna_use)) * CWCS / sumup
  num = floor(length(pro_matrix[pro_matrix != 0]) * thres)
  detail_output = matrix(data = 0, ncol = 4, nrow = num)
  sort_one_rank = order(pro_matrix, decreasing = TRUE)[1 : num]
  sort_one_prob = sort(pro_matrix, decreasing = TRUE)[1 : num]
  temp3 = ceiling(sort_one_rank/length(mirna_name))
  temp4 = sort_one_rank - (temp3 - 1) * length(mirna_name)

  for(i in 1 : num){
    detail_output[i, 1] = i
    detail_output[i, 2] = mrna_name[temp3[i]]
    detail_output[i, 3] = mirna_name[temp4[i]]
    detail_output[i, 4] = sort_one_prob[i]
  }
  colnames(detail_output) = c("rank", "Gene Symbol", "miRNA Name", "miRACLe Score")

  return(detail_output)
}

