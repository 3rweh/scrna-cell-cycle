get_rdata <- function(score_result, ordIndex)
{	
	# Loads the score results as a separate variable, uses bayes score due to better G2M separation
	# Arranges these scores with respect to the ordIndex
	bayes_cell <- score_result$bayes_score[ordIndex, ]
	# Stores each phase score as a new parameter
	g1_score <- bayes_cell$G1.score
	s_score <- bayes_cell$S.score
	g2m_score <- bayes_cell$G2M.score
	
	n_1 <- 5
	n_2 <- 4
	
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
	
	get_min_reg <- function(input, n_val)
	  # Returns the max segment region out of n for the input vector 
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
	  min_reg <- which.min(reg_vec)
	  return (min_reg)
	}
	
	get_g <- function(inp, v_2, v_3, reg)
	{
		if (50 >= length(ordIndex))
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
		    
		    d1 <- sqrt(sum((inp_b_num - v_2_b_num) ^ 2))
		    d2 <- 2*sqrt(sum((inp_b_num - v_3_b_num) ^ 2))
		    
		    d_tot <- d1+d2
		    d_vec <- append(d_vec, d_tot) 
		  }
		  
		  # Gets the indices of the n_1/ largest euclidian distances
		  r_fin_ind <- which(d_vec >= sort(d_vec, decreasing=T)[n_1/2], arr.ind=TRUE)
		  
		  fin_v <- c()
		  for (i in 1:length(r_fin_ind))
		  {
		    r_fin_n <- as.numeric(inp_splt[[brv_ind[[r_fin_ind[[i]]]]]])
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
		  
		  # Retrieves the optimal region for the input vector 
		  inp_b <- inp_splt[[reg]]
		  v_2_b <- v_2_splt[[reg]]
		  v_3_b <- v_3_splt[[reg]]
		  
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
			
			  d1 <- sqrt(sum((inp_b_num - v_2_b_num) ^ 2))
			  d2 <- 2*sqrt(sum((inp_b_num - v_3_b_num) ^ 2))
			  
			  d_tot <- d1+d2
			  d_vec <- append(d_vec, d_tot) 
		  }
		  
		  # Gets the indices of the n largest euclidian distances
		  brv_ind <- which(d_vec >= sort(d_vec, decreasing=T)[length(inp_b_splt)/2], arr.ind=TRUE)
		  
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
	
	get_s <- function(s_score, g1_score, g2m_score)
	{
	  # Function for finding the S-phase of the cycle
	  # inp = s-score vector
	  # v_2 = g1-score vector
	  # v_3= g2m_score vector
	  
	  g1_maxr <- get_max_reg(g1_score, n_2*2)
	  g1_minr <- get_min_reg(g1_score, n_2*2)
	  
	  # Forward direction
	  if(g1_minr > g1_maxr){
	    # Splits all vectors into n_2*4 regions
	    # Checks the latter 2/3 of the cycle, as the minimum is likely to be in this region
	    s_splt <- rev(split(s_score, cut(seq(s_score), n_2*4, labels = FALSE)))
	    g1_splt <- rev(split(g1_score, cut(seq(g1_score), n_2*4, labels = FALSE)))
	    g2m_splt <- rev(split(g2m_score, cut(seq(g2m_score), n_2*4, labels = FALSE)))
	    
	    # Calculates the euclidian distance for the g1 & g2 phase scores based on region
	    # In the s-phase, the distance for these two should be minimized
	    d_vec <- c()
	    for (i in 1:(length(s_splt)*(2/3)))
	    {
	      s_num <- as.numeric(s_splt[[i]])
	      g1_num <- as.numeric(g1_splt[[i]])
	      g2m_num <- as.numeric(g2m_splt[[i]])
	    
	      d_1 <- 3*sqrt(sum((s_num - g1_num) ^ 2))
	      d_2 <- sqrt(sum((g1_num - g2m_num) ^ 2))
	    
	      d_tot <- d_1 + d_2
        d_vec <- append(d_vec, d_tot) 
	    }
	  }
	  # Reverse
	  else if(g1_maxr > g1_minr){
	    # Splits all vectors into n_2*4 regions
	    # Checks the first 2/3 of the cycle, as thissi where thes phase should be located
	    s_splt <- split(s_score, cut(seq(s_score), n_2*4, labels = FALSE))
	    g1_splt <- split(g1_score, cut(seq(g1_score), n_2*4, labels = FALSE))
	    g2m_splt <- split(g2m_score, cut(seq(g2m_score), n_2*4, labels = FALSE))
	      
	    # Calculates the euclidian distance for the g1 & g2 phase scores based on region
	    # In the s-phase, the distance for these two should be minimized
	    d_vec <- c()
	    for (i in 1:(length(s_splt)*(2/3)))
	    {
	      s_num <- as.numeric(s_splt[[i]])
	      g1_num <- as.numeric(g1_splt[[i]])
	      g2m_num <- as.numeric(g2m_splt[[i]])
	      
	      d_1 <- 3*sqrt(sum((s_num - g1_num) ^ 2))
	      d_2 <- sqrt(sum((g1_num - g2m_num) ^ 2))
	      
	      d_tot <- d_1 + d_2
	      d_vec <- append(d_vec, d_tot)  
	    }
	  }
	  
	    # Gets the indices of the n smallest euclidian distances (ie: where g1 meets with g2 and is the closest to s)
	    brv_ind <- which.min(d_vec)
	  
	    return (as.numeric(s_splt[[brv_ind]]))
	  }
	
	g1_r <- get_max_reg(g1_score, n_1)
	s_r <- get_max_reg(s_score, n_1)
	g2m_r <- get_max_reg(g2m_score, n_1)
	
	g1_b <- get_g(g1_score, s_score, g2m_score, g1_r)
	s_b <- get_s(s_score, g1_score, g2m_score)
	g2m_b <- get_g(g2m_score, s_score, g1_score, g2m_r)
	
	# Matches region-vector against the original vector
	g1_b_r <- match(g1_b, g1_score)
	s_b_r <- match(s_b, s_score)
	g2m_b_r <- match(g2m_b, g2m_score)
	
	# Generates vector with min and max of each parameter
	g1_rdata <- c(min(g1_b_r), max(g1_b_r))
	s_rdata <- c(min(s_b_r), max(s_b_r))
	g2m_rdata <- c(min(g2m_b_r), max(g2m_b_r))
	
	# rdata <- rbind(g1_rdata, s_rdata, g2m_rdata)
	rdata <- t(data.frame(g1_rdata , s_rdata, g2m_rdata))
	return(rdata)
}