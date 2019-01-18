get_start <- function(score_result, ordIndex)
{	
  # Loads the score results as a separate variable, uses bayes score due to better G2M separation
  # Arranges these scores with respect to the ordIndex
  bayes_cell <- score_result$bayes_score[ordIndex, ]
  # Stores each phase score as a new parameter
  g1_score <- bayes_cell$G1.score
  s_score <- bayes_cell$S.score
  g2m_score <- bayes_cell$G2M.score
  
  n_1 <- 6
  
  get_min_reg <- function(input)
    # Returns the max segment region out of n1 for the input vector 
  {
    val <- -Inf
    pos <- 0
    inp_splt <- split(input, cut(seq(input), n_1, labels = FALSE))
    reg_vec <- c()
    for (i in 1:length(inp_splt))
    {
      inp_num <- as.numeric(inp_splt[[i]])
      inp_mean <- mean(inp_num)
      reg_vec <- append(reg_vec, inp_mean)
    }
    min_reg <- which.min(reg_vec)
    return (min_reg)
  }
  
  get_start_p <- function(g1_score, s_score, g2m_score)
  {
    # Used in case of a shorter time-series
    if (50 >= length(ordIndex)){
      # Walks through the min-region, checking if the s-score is larger then that of the g1_score. We know that when the g1_score passes that of the s-score, we have reached the start of the cycle
      for ( i in 1:(length(g1_score)-8)){
        
        g1_p_s <- g1_score[i]
        s_p_s <- s_score[i]
        g2m_p_s <- g2m_score[i]
        
        g1_p_s_1 <- g1_score[i+1]
        s_p_s_1 <- s_score[i+1]
        g2m_p_s_1 <- g2m_score[i+1]
        
        g1_p_s_2 <- g1_score[i+2]
        s_p_s_2 <- s_score[i+2]
        g2m_p_s_2 <- g2m_score[i+2]
        
        g1_p_s_3 <- g1_score[i+3]
        s_p_s_3 <- s_score[i+3]
        g2m_p_s_3 <- g2m_score[i+3]
        
        g1_p_s_6 <- g1_score[i+6]
        s_p_s_6 <- s_score[i+6]
        g2m_p_s_6 <- g2m_score[i+6]
        
        g1_p_s_7 <- g1_score[i+7]
        s_p_s_7 <- s_score[i+7]
        g2m_p_s_7 <- g2m_score[i+7]
        
        g1_p_s_8 <- g1_score[i+8]
        s_p_s_8 <- s_score[i+8]
        g2m_p_s_8 <- g2m_score[i+8]
        
        g1_p_m <- mean(g1_p_s + g1_p_s_1 + g1_p_s_2 + g1_p_s_3)
        s_p_m <- mean(s_p_s + s_p_s_1 + s_p_s_2 + s_p_s_3)
        g2m_p_m <- mean(g2m_p_s + g2m_p_s_1 + g2m_p_s_2 + g2m_p_s_3)
        
        g1_p_m_2 <- mean(g1_p_s_6 + g1_p_s_7 + g1_p_s_8)
        s_p_m_2 <- mean(s_p_s_6 + s_p_s_7 + s_p_s_8)
        g2m_p_m_2 <- mean(g2m_p_s_6 + g2m_p_s_7 + g2m_p_s_8)
        
        if(g1_p_s >= s_p_s && g1_p_s >= g2m_p_s && g1_p_s_1 >= s_p_s_1 && g1_p_s_1 >= g2m_p_s_1 && g1_p_m >= s_p_m && g1_p_m >= g2m_p_m && g1_p_m_2 >= s_p_m_2 && g1_p_m_2 >= g2m_p_m_2){
          start_pos <- i
          break
        }else if(i == (length(g1_score)-8)){
          start_pos <- i+8
        }
      }
      return(g1_score[[start_pos]])
    }else{
      # Splits all vectors into n_1 regions
      g1_splt <- split(g1_score, cut(seq(g1_score), n_1, labels = FALSE))
      s_splt <- split(s_score, cut(seq(s_score), n_1, labels = FALSE))
      g2m_splt <- split(g2m_score, cut(seq(g2m_score), n_1, labels = FALSE))
    
      g2m_minr <- get_min_reg(g2m_score)
      
      g1_vec <- c()
      s_vec <- c()
      g2m_vec <- c()
      for(i in 1:g2m_minr){
        g1_num <- as.numeric(g1_splt[[i]])
        s_num <- as.numeric(s_splt[[i]])
        g2m_num <- as.numeric(g2m_splt[[i]])
        
        g1_vec <- append(g1_vec, g1_num) 
        s_vec <- append(s_vec, s_num) 
        g2m_vec <- append(g2m_vec, g2m_num) 
      }
    
      # Walks through the min-region, checking if the s-score is larger then that of the g1_score. We know that when the g1_score passes that of the s-score, we have reached the start of the cycle
      for ( i in 1:(length(g1_vec)-10)){
        
        g1_p_s <- g1_vec[i]
        s_p_s <- s_vec[i]
        g2m_p_s <- g2m_vec[i]
        
        g1_p_s_1 <- g1_vec[i+1]
        s_p_s_1 <- s_vec[i+1]
        g2m_p_s_1 <- g2m_vec[i+1]
       
        g1_p_s_2 <- g1_vec[i+2]
        s_p_s_2 <- s_vec[i+2]
        g2m_p_s_2 <- g2m_vec[i+2]
        
        g1_p_s_3 <- g1_vec[i+3]
        s_p_s_3 <- s_vec[i+3]
        g2m_p_s_3 <- g2m_vec[i+3]
        
        g1_p_s_4 <- g1_vec[i+4]
        s_p_s_4 <- s_vec[i+4]
        g2m_p_s_4 <- g2m_vec[i+4]
        
        g1_p_s_5 <- g1_vec[i+5]
        s_p_s_5 <- s_vec[i+5]
        g2m_p_s_5 <- g2m_vec[i+5]
        
        g1_p_s_6 <- g1_vec[i+6]
        s_p_s_6 <- s_vec[i+6]
        g2m_p_s_6 <- g2m_vec[i+6]
        
        g1_p_s_7 <- g1_vec[i+7]
        s_p_s_7 <- s_vec[i+7]
        g2m_p_s_7 <- g2m_vec[i+7]
        
        g1_p_s_8 <- g1_vec[i+8]
        s_p_s_8 <- s_vec[i+8]
        g2m_p_s_8 <- g2m_vec[i+8]
        
        g1_p_s_9 <- g1_vec[i+9]
        s_p_s_9 <- s_vec[i+9]
        g2m_p_s_9 <- g2m_vec[i+9]
        
        g1_p_s_10 <- g1_vec[i+10]
        s_p_s_10 <- s_vec[i+10]
        g2m_p_s_10 <- g2m_vec[i+10]
        
        g1_p_m <- mean(g1_p_s_1 + g1_p_s_2 + g1_p_s_3 + g1_p_s_4)
        s_p_m <- mean(s_p_s_1 + s_p_s_2 + s_p_s_3 + s_p_s_4)
        g2m_p_m <- mean(g2m_p_s_1 + g2m_p_s_2 + g2m_p_s_3 + g2m_p_s_4)
        
        g1_p_m_2 <- mean(g1_p_s_5 + g1_p_s_6 + g1_p_s_7)
        s_p_m_2 <- mean(s_p_s_5 + s_p_s_6 + s_p_s_7)
        g2m_p_m_2 <- mean(g2m_p_s_5 + g2m_p_s_6 + g2m_p_s_7)
        
        g1_p_m_3 <- mean(g1_p_s_8 + g1_p_s_9 + g1_p_s_10)
        s_p_m_3 <- mean(s_p_s_8 + s_p_s_9 + s_p_s_10)
        g2m_p_m_3 <- mean(g2m_p_s_8 + g2m_p_s_9 + g2m_p_s_10)
        
        if(g1_p_s >= s_p_s && g1_p_s >= g2m_p_s && g1_p_m >= s_p_m && g1_p_m >= g2m_p_m && g1_p_m_2 >= s_p_m_2 && g1_p_m_2 >= g2m_p_m_2 && g1_p_m_3 >= s_p_m_3 && g1_p_m_3 >= g2m_p_m_3){
          start_pos <- i
          break
        }else if(i == (length(g1_vec)-10)){
          start_pos <- i+10
        }
      }
      return(g1_vec[[start_pos]])
    }
  }

  start_p <- get_start_p(g1_score, s_score, g2m_score)
  
  # Matches region-vector against the original vector
  start = match(start_p, g1_score)
  
  return(start)
}