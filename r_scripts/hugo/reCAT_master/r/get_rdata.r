get_rdata <- function(score_result, ordIndex)
{	
	# Loads the score results as a separate variable, uses bayes score due to better G2M separation
	# Arranges these scores with respect to the ordIndex
	bayes_cell <- score_result$bayes_score[ordIndex, ]
	# Stores each phase score as a new parameter
	g1_score <- bayes_cell$G1.score
	s_score <- bayes_cell$S.score
	g2m_score <- bayes_cell$G2M.score
	
	get_max_reg <- function(input)
	# Returns the max segment region out of 6 for the input vector 
	{
		val = -Inf
		pos = 0
		inp_splt <- split(input, cut(seq(input), 6, labels = FALSE))
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
	
	get_sep <- function(inp, v_2, v_3, reg)
	{
		val = Inf
		pos = 0
		n_1 = 6
		n_2 = 6
		
		
		# Splits all vectors into fours region
		inp_splt <- split(inp, cut(seq(inp), n_1, labels = FALSE))
		v_2_splt <- split(v_2, cut(seq(v_2), n_1, labels = FALSE))
		v_3_splt <- split(v_3, cut(seq(v_3), n_1, labels = FALSE))
		
		# Retrieves the optimal region for the input vector 
		inp_b = inp_splt[[reg]]
		v_2_b = v_2_splt[[reg]]
		v_3_b = v_3_splt[[reg]]
		
		# Splits the best region into 6 segments
		inp_b_splt <- split(inp_b, cut(seq(inp_b), n_2, labels = FALSE))
		v_2_b_splt <- split(v_2_b, cut(seq(v_2_b), n_2, labels = FALSE))
		v_3_b_splt <- split(v_3_b, cut(seq(v_3_b), n_2, labels = FALSE)) 
		
		
		
		d_vec <- c()
		for (i in 1:length(inp_b_splt))
		{
			# Calculates euclidian distance for each section
			inp_b_num = as.numeric(inp_b_splt[[i]])
			v_2_b_num = as.numeric(v_2_b_splt[[i]])
			v_3_b_num = as.numeric(v_3_b_splt[[i]])
			
			d1 = sqrt(sum((inp_b_num - v_2_b_num) ^ 2))
			d2 = sqrt(sum((inp_b_num - v_3_b_num) ^ 2))
			d_tot = d1+d2
			d_vec <- append(d_vec, d_tot) 
		}
		max_reg <- which.max(d_vec)
		return (as.numeric(inp_b_splt[[max_reg]]))
	}
	
	g1_r <- get_max_reg(g1_score)
	s_r <- get_max_reg(s_score)
	g2m_r <- get_max_reg(g2m_score)
	
	g1_b <- get_sep(g1_score, s_score, g2m_score, g1_r)
	s_b <- get_sep(s_score, g1_score, g2m_score, s_r)
	g2m_b <- get_sep(g2m_score, g1_score, s_score, g2m_r)
	
	# Matches region-vector against the original vector
	g1_b_r = match(g1_b, g1_score)
	s_b_r = match(s_b, s_score)
	g2m_b_r = match(g2m_b, g2m_score)
	
	# Generates vector with min and max of each parameter
	g1_rdata <- c(min(g1_b_r), max(g1_b_r))
	s_rdata <- c(min(s_b_r), max(s_b_r))
	g2m_rdata = c(min(g2m_b_r), max(g2m_b_r))
	
	# rdata <- rbind(g1_rdata, s_rdata, g2m_rdata)
	rdata = t(data.frame(g1_rdata , s_rdata, g2m_rdata))
	return(rdata)
}