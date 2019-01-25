# Made by: Hugo Swenson, 2019-01-25

get_rdata <- function(score_result, ordIndex)
{	
	# Loads the score results as a separate variable, uses bayes score due to better G2M separation
	# Arranges these scores with respect to the ordIndex
	bayes_cell <- score_result$bayes_score[ordIndex, ]
	# Stores each phase score as a new parameter
	g1_score <- bayes_cell$G1.score
	s_score <- bayes_cell$S.score
	g2m_score <- bayes_cell$G2M.score
	
	n_1 = 4
	
	get_max_reg <- function(input, n_val)
	# Returns the max segment region out of n1 for the input vector 
	{
		val <- -Inf
		pos <- 0
		inp_splt <- split(input, cut(seq(input), n_val, labels = FALSE))
		reg_vec <- c()
		for (i in 1:length(inp_splt))
		{
			inp_num <- as.numeric(inp_splt[[i]])
			inp_mean <- mean(inp_num)
			reg_vec <- append(reg_vec, inp_mean)
		}
		max_reg <- which.max(reg_vec)
		return (max_reg)
	}
	
	get_g <- function(inp, v_2, v_3, reg)
	{
		if (100 >= length(ordIndex))
		{
		  # Splits all vectors into n_1 regions
		  inp_splt <- split(inp, cut(seq(inp), n_1, labels = FALSE))
		  v_2_splt <- split(v_2, cut(seq(v_2), n_1, labels = FALSE))
		  v_3_splt <- split(v_3, cut(seq(v_3), n_1, labels = FALSE))
		  
		  # Calculates the mean value of the n_1 regions
		  br_vec <- c()
		  for (i in 1:length(inp_splt))
		  {
		    inp_nbr <- as.numeric(inp_splt[[i]])
		    inp_mbr <- mean(inp_nbr)
		    br_vec <- append(br_vec, inp_mbr)
		  }
		  
		  # Gets the indices of the n_1/2 largest mean
		  brv_ind <- which(br_vec >= sort(br_vec, decreasing=T)[n_1/2], arr.ind=TRUE)
		  
		  # Calculates euclidian distance for each of the max value mean regions
		  d_vec <- c()
		  for (i in 1:length(brv_ind))
		  {
		    inp_b_num <- as.numeric(inp_splt[[(brv_ind[[i]])]])
		    v_2_b_num <- as.numeric(v_2_splt[[(brv_ind[[i]])]])
		    v_3_b_num <- as.numeric(v_3_splt[[(brv_ind[[i]])]])
		    
		    inp_b_num[inp_b_num < max(inp_b_num)]
		    v_2_b_num[v_2_b_num < max(v_2_b_num)]
		    v_2_b_num[v_3_b_num < max(v_2_b_num)]
		    
		    d1 <- sqrt(sum((inp_b_num - v_2_b_num) ^ 2))
		    d2 <- sqrt(sum((inp_b_num - v_3_b_num) ^ 2))
		    
		    # Control sequence so that the distance is not due to the g1-phase separating from the g2m-phase
		    if (mean(v_3_b_num) > mean(inp_b_num)){
		      d2 <- -1*d2
		    }
		    
		    d_tot <- d1+d2
		    d_vec <- append(d_vec, d_tot) 
		  }
		  
		  # Gets the indices of the n_1/ largest euclidian distances
		  r_fin_ind <- which(d_vec >= sort(d_vec, decreasing=T)[n_1/2], arr.ind=TRUE)
		  
		  fin_v <- c()
		  for (i in 1:length(r_fin_ind))
		  {
		    r_fin_n <- as.numeric(inp_splt[[brv_ind[[r_fin_ind[[i]]]]]])
		    r_fin_n[r_fin_n < max( r_fin_n)]
		    r_fin_m <- mean(r_fin_n)
		    fin_v <- append(fin_v, r_fin_m)
		  }
		  
		  # Returns the index no. of the largest mean value of the n max euclidian distances
		  max_reg <- which.max(fin_v)
		  
		  return (as.numeric(inp_splt[[(brv_ind[[r_fin_ind[[max_reg]]]])]]))
		} else {
		  
		  # Splits all vectors into n_1 regions
		  inp_splt <- split(inp, cut(seq(inp), n_1, labels = FALSE))
		  v_2_splt <- split(v_2, cut(seq(v_2), n_1, labels = FALSE))
		  v_3_splt <- split(v_3, cut(seq(v_3), n_1, labels = FALSE))
		  
		  # Gets the optimal region
		  inp_b <- inp_splt[[reg]]
		  v_2_b <- v_2_splt[[reg]]
		  v_3_b <- v_3_splt[[reg]]
		  
		  # Obtains the closest integer 8rounded upward) for the length of the score vector divided by 21
		  n_2 = ceiling(length(inp_b)/21)

		  # Splits the best region into n_2 segments
		  inp_b_splt <- split(inp_b, cut(seq(inp_b), n_2, labels = FALSE))
		  v_2_b_splt <- split(v_2_b, cut(seq(v_2_b), n_2, labels = FALSE))
		  v_3_b_splt <- split(v_3_b, cut(seq(v_3_b), n_2, labels = FALSE))
		  
		  # Calculates euclidian distance for each region
		  d_vec <- c()
		  for (i in 1:length(inp_b_splt))
		  {
			  inp_b_num <- as.numeric(inp_b_splt[[i]])
			  v_2_b_num <- as.numeric(v_2_b_splt[[i]])
			  v_3_b_num <- as.numeric(v_3_b_splt[[i]])
			  
			  d1 <- 0.5*sqrt(sum((inp_b_num - v_2_b_num) ^ 2))
			  d2 <- sqrt(sum((inp_b_num - v_3_b_num) ^ 2))
			  
			  # Control sequence so that the distance is not due to the g1-phase separating from the g2m-phase
			  if (mean(v_3_b_num) > mean(inp_b_num)){
			    d2 <- -1*d2
			  }
			
			  d_tot <- d1+d2
			  d_vec <- append(d_vec, d_tot) 
		  }
		  
		  # Gets the indices of the n largest euclidian distances
		  brv_ind <- which(d_vec >= sort(d_vec, decreasing=T)[ceiling(length(inp_b_splt)/2)], arr.ind=TRUE)
		  
		  # Calculates the mean value of the n regions
		  br_vec <- c()
		  for (i in 1:length(brv_ind))
		  {
		    inp_nbr <- as.numeric(inp_b_splt[[(brv_ind[[i]])]])
		    inp_mbr <- mean(inp_nbr)
		    br_vec <- append(br_vec, inp_mbr)
		  }
		  
		  # Returns the largest mean value of the n max euclidian distances
		  max_reg <- which.max(br_vec)
		  return (as.numeric(inp_b_splt[[(brv_ind[[max_reg]])]]))
		}
	}
	
	get_s <- function(s_score, g1_score, g2m_score, g1_s)
	{
	  if (50 >= length(ordIndex))
	  {
	    n_2 <- ceiling(length(s_score)/4)
	    
	    s_splt <- split(s_score, cut(seq(s_score), n_2, labels = FALSE))
	    g1_splt <- split(g1_score, cut(seq(g1_score), n_2, labels = FALSE))
	    g2m_splt <- split(g2m_score, cut(seq(g2m_score), n_2, labels = FALSE))
	    
	    # Calculates the euclidian distance for the g1 & g2 phase scores based on region
	    # In the s-phase, the distance for these two should be minimized
	    d_vec <- c()
	    for (i in 1:(length(s_splt)))
	    {
	      d_tot = 0
	      s_num <- as.numeric(s_splt[[i]])
	      g1_num <- as.numeric(g1_splt[[i]])
	      g2m_num <- as.numeric(g2m_splt[[i]])
	      
	      d_1 <- sqrt(sum((s_num - g1_num) ^ 2))
	      d_2 <- sqrt(sum((s_num - g2m_num) ^ 2))
	      d_3 <- sqrt(sum((g1_num - g2m_num) ^ 2)) 
	      
	      # Control sequence so that the distance is not due to the g1-phase separating from the g2m-phase
	      if (mean(g2m_num) > mean(s_num)){
	        d_2 <- 2*d_2
	      } else if (mean(g1_num) > mean(s_num)){
	        d_1 <- 2*d_1
	      }
	      
	      d_tot <- d_1 + d_2 + d_3
	      d_vec <- append(d_vec, d_tot) 
	    }
	    # Gets the indices of the n smallest euclidian distances (ie: where g1 meets with g2 and is the closest to s)
	    brv_ind <- which.min(d_vec)
	    return (as.numeric(s_splt[[brv_ind]]))
	  }else{
	    
	    # Creates a new vecotr ranging from the min of the best region for the G1 phase, to the end of the vector
	    new_s <- s_score[g1_s:length(ordIndex)]
	    new_g1 <- g1_score[g1_s:length(ordIndex)]
	    new_g2m <- g2m_score[g1_s:length(ordIndex)]
	    
	    # Obtains the closest integer value rounded upwards, for the length of the vector divided by 21
	    n_2 <- ceiling(length(new_s)/21)
	  
	    s_splt <- split(new_s, cut(seq(new_s), n_2, labels = FALSE))
	    g1_splt <- split(new_g1, cut(seq(new_g1), n_2, labels = FALSE))
	    g2m_splt <- split(new_g2m, cut(seq(new_g2m), n_2, labels = FALSE))
	    
	    # Calculates the euclidian distance for the g1 & g2 phase scores based on region
	    # In the s-phase, the distance for these two should be minimized
	    d_vec <- c()
	    for (i in 1:(length(s_splt)))
	    {
	      d_tot = 0
	      s_num <- as.numeric(s_splt[[i]])
	      g1_num <- as.numeric(g1_splt[[i]])
	      g2m_num <- as.numeric(g2m_splt[[i]])
	    
	      d_1 <- 3*sqrt(sum((s_num - g1_num) ^ 2))
	      d_2 <- sqrt(sum((s_num - g2m_num) ^ 2))
	      d_3 <- sqrt(sum((g1_num - g2m_num) ^ 2))
	      d_tot <- d_1 + d_2 + d_3
        d_vec <- append(d_vec, d_tot) 
	    }
	    # Gets the indices of the n smallest euclidian distances (ie: where g1 meets with g2 and is the closest to s)
	    brv_ind <- which.min(d_vec)
	    return (as.numeric(s_splt[[brv_ind]]))
	  }
	}
	
	# Obtains the best regions for each phase
	g1_r <- get_max_reg(g1_score, n_1)
	s_r <- get_max_reg(s_score, n_1)
	g2m_r <- get_max_reg(g2m_score, n_1)
	
	# Obtains the phase regions for g1 and g2m
	g1_b <- get_g(g1_score, s_score, g2m_score, g1_r)
	g2m_b <- get_g(g2m_score, s_score, g1_score, g2m_r)
	
	# Matches region-vector against the original vector
	g1_b_r <- match(g1_b, g1_score)
	g2m_b_r <- match(g2m_b, g2m_score)
	
	# Obtains the s-phase region
	g1_s <- min(g1_b_r)
	s_b <- get_s(s_score, g1_score, g2m_score, g1_s)
	s_b_r <- match(s_b, s_score)
	
	# Generates vector with min and max of each parameter
	g1_rdata <- c(min(g1_b_r), max(g1_b_r))
	s_rdata <- c(min(s_b_r), max(s_b_r))
	g2m_rdata <- c(min(g2m_b_r), max(g2m_b_r))
	
	# rdata <- rbind(g1_rdata, s_rdata, g2m_rdata)
	rdata <- t(data.frame(g1_rdata , s_rdata, g2m_rdata))
	return(rdata)
}