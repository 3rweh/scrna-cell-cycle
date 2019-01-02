get_start <- function(score_result, ordIndex)
{	
  # Loads the score results as a separate variable, uses bayes score due to better G2M separation
  # Arranges these scores with respect to the ordIndex
  bayes_cell <- score_result$bayes_score[ordIndex, ]
  # Stores each phase score as a new parameter
  g1_score <- bayes_cell$G1.score
  s_score <- bayes_cell$S.score
  g2m_score <- bayes_cell$G2M.score
  
  n_1 = 6
  
  get_min_reg <- function(input)
    # Returns the max segment region out of n1 for the input vector 
  {
    val = -Inf
    pos = 0
    inp_splt <- split(input, cut(seq(input), n_1, labels = FALSE))
    reg_vec <- c()
    for (i in 1:length(inp_splt))
    {
      inp_num <- as.numeric(inp_splt[[i]])
      inp_mean = mean(inp_num)
      reg_vec <- append(reg_vec, inp_mean)
    }
    min_reg <- which.min(reg_vec)
    return (min_reg)
  }
  
  get_max_reg <- function(input)
    # Returns the max segment region out of n1 for the input vector 
  {
    val = -Inf
    pos = 0
    inp_splt <- split(input, cut(seq(input), n_1, labels = FALSE))
    reg_vec <- c()
    for (i in 1:length(inp_splt))
    {
      inp_num <- as.numeric(inp_splt[[i]])
      inp_mean = mean(inp_num)
      reg_vec <- append(reg_vec, inp_mean)
    }
    max_reg <- which.max(reg_vec)
    return (max_reg)
  }
  
  get_start_p <- function(g1_score, s_score, g2m_score)
  {
    # Gets the region where g2m is the smallest (ie, the peak of the g1 phase)
    g2m_minr <- get_min_reg(g2m_score)
    # Gets the region where g2m is the largest (ie, the min of the g1 phase)
    g2m_maxr <- get_max_reg(g2m_score)
  
  
    # Checks the direction of the cycle, if the peak of the g2m phase is in front of the min, it is forward
    if(g2m_maxr > g2m_minr){
      # Splits all vectors into n_1 regions
      g1_splt <- split(g1_score, cut(seq(g1_score), n_1, labels = FALSE))
      s_splt <- split(s_score, cut(seq(s_score), n_1, labels = FALSE))
      g2m_splt <- split(g2m_score, cut(seq(g2m_score), n_1, labels = FALSE))
    
      s_vec <- c()
      for(i in 1:g2m_minr){
        g1_num = as.numeric(g1_splt[[i]])
        s_num = as.numeric(s_splt[[i]])
        g2m_num = as.numeric(g2m_splt[[i]])
        
        d_1 = sqrt(sum((g1_num - s_num) ^ 2))
        d_2 = sqrt(sum((g1_num - g2m_num) ^ 2))
        d_3 = sqrt(sum((s_num - g2m_num) ^ 2)) 
        d_tot = d_1 + d_2 + d_3
        s_vec <- append(s_vec, d_tot) 
      }
      # Gets the indices of the n smallest euclidian distances (ie: where g1 meets with g2 and is the closest to s)
      start_ind <- which.min(s_vec)
    
      g1_s_reg = as.numeric(g1_splt[[start_ind]])
      s_s_reg = as.numeric(s_splt[[start_ind]])
    
      # Walks through the min-region, checking if the s-score is larger then that of the g1_score. We know that when the g1_score passes that of the s-score, we have reached the start of the cycle
      for ( i in 1:length(g1_s_reg)){
        g1_p_s = g1_s_reg[i]
        s_p_s = s_s_reg[i]
        if(g1_p_s > s_p_s){
          start_pos = i
          break
        } else {
          # If the G1 phase score is always higher then that of the S-phase score,  that means that the start is located in position 1
          start_pos = 1 
        }
      }
      return(g1_s_reg[[start_pos]])
    }else if(g2m_minr > g2m_maxr){
      # Splits all vectors into n_1 regions
      g1_splt <- rev(split(g1_score, cut(seq(g1_score), n_1, labels = FALSE)))
      s_splt <- rev(split(s_score, cut(seq(s_score), n_1, labels = FALSE)))
      g2m_splt <- rev(split(g2m_score, cut(seq(g2m_score), n_1, labels = FALSE)))
    
      s_vec <- c()
      for(i in length(g1_splt):g2m_minr){
        g1_num = as.numeric(g1_splt[[i]])
        s_num = as.numeric(s_splt[[i]])
        g2m_num = as.numeric(g2m_splt[[i]])
      
        d_1 = sqrt(sum((g1_num - s_num) ^ 2))
        d_2 = sqrt(sum((g1_num - g2m_num) ^ 2))
        d_3 = sqrt(sum((s_num - g2m_num) ^ 2)) 
        d_tot = d_1 + d_2 + d_3
        s_vec <- append(s_vec, d_tot) 
      }
      # Gets the indices of the n smallest euclidian distances (ie: where g1 meets with g2 and is the closest to s)
      start_ind <- which.min(s_vec)
    
      g1_s_reg = as.numeric(g1_splt[[start_ind]])
      s_s_reg = as.numeric(s_splt[[start_ind]])
    
      # Walks through the min-region, checking if the s-score is larger then that of the g1_score. We know that when the g1_score passes that of the s-score, we have reached the start of the cycle
      for (i in length(g1_s_reg):1){
        g1_p_s = g1_s_reg[i]
        s_p_s = s_s_reg[i]
        if(g1_p_s > s_p_s){
          start_pos = i
          break
        }else{
          # If the G1 phase score is always higher then that of the S-phase score,  that means that the start is located in position 1
          start_pos = length(g1_s_reg)
        }
      }
    }
    return(g1_s_reg[[start_pos]])
  }

  start_p <- get_start_p(g1_score, s_score, g2m_score)
  
  # Matches region-vector against the original vector
  start = match(start_p, g1_score)
  
  return(start)
}