#' @title miracle is the main function
#'
#' @param seqScore 
#' @param mirtarList 
#' @param mirExpr 
#' @param tarExpr 
#' @param exprFilter 
#' @param OutputSelect 
#' @param ID_pos 
#'
#' @return
#' @export
#'
#' @examples
miracle = function(seqScore, mirtarList, mirExpr, tarExpr, exprFilter = 1, OutputSelect = TRUE, ID_pos = c()){
  input_list = matrix_transfer(seqScore)
  x = name_func(mirtarList, mirExpr, tarExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(tarExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], exprFilter)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  if(length(ID_pos) == 0){
    z9 = miracle_detail(z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], OutputSelect = TRUE, select_ID = mirtarList[, 1])
  }
  else{
    z9 = miracle_detail(z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], OutputSelect = TRUE, select_ID = ID_pos)
  }
  
  return(z9)
}


#' @title miracle_mirna is used to get the ranked result of the targets of a given mirna
#'
#' @param target_name 
#' @param seqScore 
#' @param mirtarList 
#' @param mirExpr 
#' @param tarExpr 
#' @param exprFilter 
#' @param ID_pos 
#'
#' @return
#' @export
#'
#' @examples
miracle_mirna = function(target_name, seqScore, mirtarList, mirExpr, tarExpr, exprFilter = 1, ID_pos = c()){
  input_list = matrix_transfer(seqScore)
  x = name_func(mirtarList, mirExpr, tarExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(tarExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], exprFilter)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  if(length(ID_pos) == 0){
    z9 = mirna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = mirtarList[, 1])
  }
  else{
    z9 = mirna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = ID_pos)
  }
  return(z9)
}



#' @title miracle_target is used to get the ranked result of the targets of a given mrna
#'
#' @param target_name 
#' @param seqScore 
#' @param mirtarList 
#' @param mirExpr 
#' @param tarExpr 
#' @param exprFilter 
#' @param ID_pos 
#'
#' @return
#' @export
#'
#' @examples
miracle_target = function(target_name, seqScore, mirtarList, mirExpr, tarExpr, exprFilter = 1, ID_pos = c()){
  input_list = matrix_transfer(seqScore)
  x = name_func(mirtarList, mirExpr, tarExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(tarExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], exprFilter)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  if(length(ID_pos) == 0){
    z9 = mrna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = mirtarList[, 1])
  }
  else{
    z9 = mrna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = ID_pos)
  }
  return(z9)
}


#' @title miracle_multi_mirna is used to get the result of the prediction result of multiple given mirnas
#'
#' @param target_name 
#' @param seqScore 
#' @param mirtarList 
#' @param mirExpr 
#' @param tarExpr 
#' @param exprFilter 
#' @param ID_pos 
#'
#' @return
#' @export
#'
#' @examples
miracle_multi_mirna = function(target_name, seqScore, mirtarList, mirExpr, tarExpr, exprFilter = 1, ID_pos = c()){
  input_list = matrix_transfer(seqScore)
  x = name_func(mirtarList, mirExpr, tarExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(tarExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], exprFilter)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  if(length(ID_pos) == 0){
    z9 = multi_mirna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = mirtarList[, 1])
  }
  else{
    z9 = multi_mirna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = ID_pos)
  }
  return(z9)
}

#' @title miracle_multi_target is used to get the result of the prediction result of multiple given mrnas
#'
#' @param target_name 
#' @param seqScore 
#' @param mirtarList 
#' @param mirExpr 
#' @param tarExpr 
#' @param exprFilter 
#' @param ID_pos 
#'
#' @return
#' @export
#'
#' @examples
miracle_multi_target = function(target_name, seqScore, mirtarList, mirExpr, tarExpr, exprFilter = 1, ID_pos = c()){
  input_list = matrix_transfer(seqScore)
  x = name_func(mirtarList, mirExpr, tarExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(tarExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], exprFilter)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  if(length(ID_pos) == 0){
    z9 = multi_mrna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = mirtarList[, 1])
  }
  else{
    z9 = multi_mrna_extract(target_name, z3[[1]], z3[[2]], z7[[3]], z8, mirtarList[, 1], z3[[3]], z3[[4]], select_ID = ID_pos)
  }
  return(z9)
}





#matrix_transfer is used to get the sequence matching matrix
matrix_transfer = function(input_matrix){
  mirna_name = unique(input_matrix[, 2])
  mrna_name = unique(input_matrix[, 1])
  mod_matrix = matrix(data = 0, ncol = length(mirna_name), nrow = length(mrna_name))
  pos1 = match(input_matrix[, 2], mirna_name)
  pos2 = match(input_matrix[, 1], mrna_name)
  for(i in 1:nrow(input_matrix)){
    mod_matrix[pos2[i], pos1[i]] = as.numeric(input_matrix[i, 3])
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
  mirna_use = mirna_base[mirna_name1%in%mirna, mirna_name]
  mirna_name1 = mirna_name1[mirna_name1%in%mirna]
  return(list(mirna_use, mirna_name1))
}
mrna_matrix = function(mrna_base, mrna_name, mrna){
  mrna_exist = rep(0, nrow(mrna_base))
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist, mrna_name]
  mrna_name2 = mrna_name1[mrna_exist]
  return(list(mrna_use, mrna_name2))
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

#the following function sum_miracle and miracle_detail are the core functions of the miracle algorithm.
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
#The default is the first percent pairs, while the users can set the OutputSelect to be False to get full result
miracle_detail = function(mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, OutputSelect = TRUE, select_ID){
  if(OutputSelect == TRUE){
    thres = 0.01
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
    mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
    mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
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
#mirna_extract is used to get the ranked result of a specific mirna
mirna_extract = function(target_name, mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, select_ID){
  ID_poss = match(select_ID, ID)
  ID_pos = ID_poss[which(!is.na(ID_poss))]
  if(length(ID_pos) == 0){
    print("The selected samples are not in the predicted list.")
    return("The selected samples are not in the predicted list.")
  }
  else{
    if(target_name%in%mirna_name == FALSE){
      print("This mirna is not in the predicted list.")
      return("This mirna is not in the predicted list.")
    }
    else{
      mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
      mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
      mirna_use = mirna_use[, ID_pos]
      mrna_use = mrna_use[, ID_pos]
      ID_new = ID[ID_pos]
      pos = which(mirna_name == target_name)
      num = length(which(CWCS[pos, ] != 0))
      if(is.vector(mirna_use) == TRUE){
        mirna_use = as.matrix(mirna_use)
      }
      if(is.vector(mrna_use) == TRUE){
        mrna_use = as.matrix(mrna_use)
      }

      mrna_ind = matrix(data = 0, ncol = (3 * ncol(mrna_use)), nrow = num, byrow = TRUE)
      nonzero_count = c()
      mrna_ind_name = c()
      pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
      rank_num = seq(1, num, 1)
      for(m in 1:ncol(mirna_use)){
        mrna_ind_name = c(mrna_ind_name, ID_new[m], "Gene Symbol", "miRACLe Score")
        pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
        pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
        mrna_list = pro_matrix_temp[pos,]
        nonzero_count = c(nonzero_count, length(which(mrna_list != 0)))
        rank = order(mrna_list, decreasing = TRUE)
        mrna_ind[ , (3 * m - 2)] = rank_num
        mrna_ind[ , (3 * m - 1)] = mrna_name[rank[1:num]]
        mrna_ind[ , (3 * m)] = mrna_list[rank[1:num]]
      }
      colnames(mrna_ind) = mrna_ind_name
      mrna_ind = mrna_ind[1:max(nonzero_count), ]
      pro_matrix_agg = pro_matrix_agg/ncol(mirna_use)
      mrna_list_agg = pro_matrix_agg[pos, ]
      rank = order(mrna_list_agg, decreasing = TRUE)
      mrna_name_group = mrna_name[rank[1:num]]
      score_group = mrna_list_agg[rank[1:num]]
      mrna_group = cbind(rank_num, mrna_name_group, score_group)
      colnames(mrna_group) = c("rank", "Gene Symbol", "miRACLe Score")
      mrna_group = mrna_group[1:length(which(mrna_list_agg != 0)),]
    }
  }
  return(list(mrna_ind, mrna_group))
}
#mrna_extract is used to get the ranked result of a specific mrna
mrna_extract = function(target_name, mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, select_ID){
  ID_poss = match(select_ID, ID)
  ID_pos = ID_poss[which(!is.na(ID_poss))]
  if(length(ID_pos) == 0){
    print("The selected samples are not in the predicted list.")
    return("The selected samples are not in the predicted list.")
  }
  else{
    if(target_name%in%mrna_name == FALSE){
      print("This mrna is not in the predicted list.")
      return("This mrna is not in the predicted list.")
    }
    else{
      mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
      mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
      mirna_use = mirna_use[, ID_pos]
      mrna_use = mrna_use[, ID_pos]
      ID_new = ID[ID_pos]
      pos = which(mrna_name == target_name)
      num = length(which(CWCS[, pos] != 0))
      if(is.vector(mirna_use) == TRUE){
        mirna_use = as.matrix(mirna_use)
      }
      if(is.vector(mrna_use) == TRUE){
        mrna_use = as.matrix(mrna_use)
      }
      mirna_ind = matrix(data = 0, ncol = (3 * ncol(mirna_use)), nrow = num, byrow = TRUE)
      rank_num = seq(1, num, 1)
      nonzero_count = c()
      mirna_ind_name = c()
      pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
      for(m in 1:ncol(mirna_use)){
        mirna_ind_name = c(mirna_ind_name, ID_new[m], "miRNA Name", "miRACLe Score")
        pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
        pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
        mirna_list = pro_matrix_temp[, pos]
        nonzero_count = c(nonzero_count, length(which(mirna_list != 0)))
        rank = order(mirna_list, decreasing = TRUE)
        mirna_ind[ , (3 * m - 2)] = rank_num
        mirna_ind[ , (3 * m - 1)] = mirna_name[rank[1:num]]
        mirna_ind[ , (3 * m)] = mirna_list[rank[1:num]]
      }
      colnames(mirna_ind) = mirna_ind_name
      mirna_ind = mirna_ind[1:max(nonzero_count), ]
      mirna_list_agg = pro_matrix_agg[, pos]
      rank = order(mirna_list_agg, decreasing = TRUE)
      mirna_name_group = mirna_name[rank[1:num]]
      score_group = mirna_list_agg[rank[1:num]]
      mirna_group = cbind(rank_num, mirna_name_group, score_group)
      colnames(mirna_group) = c("rank", "miRNA Name", "miRACLe Score")
      mirna_group = mirna_group[1:length(which(mirna_list_agg != 0)), ]
    }
  }
  return(list(mirna_ind, mirna_group))
}
#multi_mirna_extract is used to get the pair information for multiple interested miRNAs
multi_mirna_extract = function(target_name, mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, select_ID){
  ID_poss = match(select_ID, ID)
  ID_pos = ID_poss[which(!is.na(ID_poss))]
  if(length(ID_pos) == 0){
    print("The selected samples are not in the predicted list.")
    return("The selected samples are not in the predicted list.")
  }
  else{
    mirna_poss = match(target_name, mirna_name)
    mirna_pos = mirna_poss[which(!is.na(mirna_poss))]
    if(length(mirna_pos) == 0){
      print("The selected miRNAs are not in the predicted list.")
      return("The selected miRNAs are not in the predicted list.")
    }
    else{
      mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
      mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
      mirna_use = mirna_use[, ID_pos]
      mrna_use = mrna_use[, ID_pos]
      mirna_name_select = mirna_name[mirna_pos]
      if(is.vector(mirna_use) == TRUE){
        mirna_use = as.matrix(mirna_use)
      }
      if(is.vector(mrna_use) == TRUE){
        mrna_use = as.matrix(mrna_use)
      }
      numm = length(mirna_pos) * nrow(mrna_use)
      #pro_matrix_agg saves the probability matrix of the integrated result
      pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
      #detail_output saves the final result
      detail_output = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = numm, byrow = TRUE)
      temp_name = c()
      num_max = 0
      
      for(m in 1:ncol(mirna_use)){
        pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
        pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
        pro_select = pro_matrix_temp[mirna_pos,]
        num1 = length(pro_select[pro_select != 0])
        num_max = max(num_max, num1)
        sort_rank = order(pro_select, decreasing = TRUE)[1 : num1]
        sort_prob = sort(pro_select, decreasing = TRUE)[1 : num1]
        temp3 = ceiling(sort_rank/length(mirna_pos))
        temp4 = sort_rank - (temp3 - 1) * length(mirna_pos)
        for(i in 1 : num1){
          detail_output[i, (4 * m - 3)] = i
          detail_output[i, (4 * m - 2)] = mrna_name[temp3[i]]
          detail_output[i, (4 * m - 1)] = mirna_name_select[temp4[i]]
          detail_output[i, (4 * m)] = sort_prob[i]
        }
        temp_name = c(temp_name, ID[ID_pos[m]], "Gene Symbol", "miRNA Name", "miRACLe Score")
      }
      colnames(detail_output) = temp_name
      detail_output1 = detail_output[1:num_max, ]
      
      pro_matrix_agg = pro_matrix_agg/ncol(mirna_use)
      pro_agg = pro_matrix_agg[mirna_pos,]
      num2 = length(pro_agg[pro_agg != 0])
      out_matrix = matrix(data = 0, ncol = 4, nrow = num2, byrow = TRUE)
      temp_rank = order(pro_agg, decreasing = TRUE)[1 : num2]
      temp_prob = sort(pro_agg, decreasing = TRUE)[1 : num2]
      temp5 = ceiling(temp_rank/length(mirna_pos))
      temp6 = temp_rank - (temp5 - 1) * length(mirna_pos)
      
      for(i in 1:nrow(out_matrix)){
        out_matrix[i, 1] = i
        out_matrix[i, 2] = mrna_name[temp5[i]]
        out_matrix[i, 3] = mirna_name_select[temp6[i]]
        out_matrix[i, 4] = temp_prob[i]
      }
      colnames(out_matrix) = c("rank", "Gene Symbol", "miRNA Name", "miRACLe Score")
      return(list(detail_output1, out_matrix))
    }
  }
}
#multi_mrna_extract is used to get the pair information for multiple interested mRNAs
multi_mrna_extract = function(target_name, mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, select_ID){
  ID_poss = match(select_ID, ID)
  ID_pos = ID_poss[which(!is.na(ID_poss))]
  if(length(ID_pos) == 0){
    print("The selected samples are not in the predicted list.")
    return("The selected samples are not in the predicted list.")
  }
  else{
    mrna_poss = match(target_name, mrna_name)
    mrna_pos = mrna_poss[which(!is.na(mrna_poss))]
    if(length(mrna_pos) == 0){
      print("The selected mRNAs are not in the predicted list.")
      return("The selected mRNAs are not in the predicted list.")
    }
    else{
      mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
      mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
      numm = length(mrna_pos) * nrow(mirna_use)
      mirna_use = mirna_use[, ID_pos]
      mrna_use = mrna_use[, ID_pos]
      mrna_name_select = mrna_name[mrna_pos]
      if(is.vector(mirna_use) == TRUE){
        mirna_use = as.matrix(mirna_use)
      }
      if(is.vector(mrna_use) == TRUE){
        mrna_use = as.matrix(mrna_use)
      }
      #pro_matrix_agg saves the probability matrix of the integrated result
      pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
      detail_output = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = numm, byrow = TRUE)
      temp_name = c()
      num_max = 0
      
      for(m in 1:ncol(mrna_use)){
        pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
        pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
        pro_select = pro_matrix_temp[, mrna_pos]
        num1 = length(pro_select[pro_select != 0])
        num_max = max(num_max, num1)
        sort_rank = order(pro_select, decreasing = TRUE)[1 : num1]
        sort_prob = sort(pro_select, decreasing = TRUE)[1 : num1]
        temp3 = ceiling(sort_rank/length(mirna_name))
        temp4 = sort_rank - (temp3 - 1) * length(mirna_name)
        for(i in 1 : num1){
          detail_output[i, (4 * m - 3)] = i
          detail_output[i, (4 * m - 2)] = mrna_name_select[temp3[i]]
          detail_output[i, (4 * m - 1)] = mirna_name[temp4[i]]
          detail_output[i, (4 * m)] = sort_prob[i]
        }
        temp_name = c(temp_name, ID[ID_pos[m]], "Gene Symbol", "miRNA Name", "miRACLe Score")
      }
      colnames(detail_output) = temp_name
      detail_output1 = detail_output[1:num_max, ]
      
      pro_matrix_agg = pro_matrix_agg/ncol(mirna_use)
      pro_agg = pro_matrix_agg[, mrna_pos]
      num2 = length(pro_agg[pro_agg != 0])
      out_matrix = matrix(data = 0, ncol = 4, nrow = num2, byrow = TRUE)
      temp_rank = order(pro_agg, decreasing = TRUE)[1 : num2]
      temp_prob = sort(pro_agg, decreasing = TRUE)[1 : num2]
      temp5 = ceiling(temp_rank/length(mirna_name))
      temp6 = temp_rank - (temp5 - 1) * length(mirna_name)
      
      for(i in 1:nrow(out_matrix)){
        out_matrix[i, 1] = i
        out_matrix[i, 2] = mrna_name_select[temp5[i]]
        out_matrix[i, 3] = mirna_name[temp6[i]]
        out_matrix[i, 4] = temp_prob[i]
      }
      colnames(out_matrix) = c("rank", "Gene Symbol", "miRNA Name", "miRACLe Score")
      return(list(detail_output1, out_matrix))
    }
  }
}



