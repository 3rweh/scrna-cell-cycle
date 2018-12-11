get_myord <- function(score_result, ordIndex)
{	
	# Loads the score results as a separate variable, uses bayes score due to better G2M separation
	# Arranges these scores with respect to the ordIndex
	bayes_cell <- score_result$bayes_score[ordIndex, ]
	# Stores each phase score as a new parameter
	g1_score <- bayes_cell$G1.score
	
	get_ord <- function(input)
	# Returns the max and min segment regions out of 20 regions for the input vector (ie: minima and maxima) 
	{
		val = -Inf
		pos = 0
		n1 = 20
		inp_splt <- split(input, cut(seq(input), n1, labels = FALSE))
		reg_vec <- c()
		for (i in 1:length(inp_splt))
		{
			inp_num <- as.numeric(inp_splt[[i]])
			inp_mean = mean(inp_num)
			reg_vec <- append(reg_vec, inp_mean)
		}
		# Obtains the region which has the highest mean-value
		max_reg <- which.max(reg_vec)
		start_reg <- as.numeric(inp_splt[[max_reg]])
		
		# Loops through the obtained region in pairs of three, and indentifies the sub-region with the highest mean-value
		start_vec <- c()
		for (i in 1:length(start_reg)-2)
		{
			inp_num_1 <- as.numeric(inp_splt[[i]])
			inp_num_2 <- as.numeric(inp_splt[[i+1]])
			inp_num_3 <- as.numeric(inp_splt[[i+2]])
			num_vec <- c(inp_num_1, inp_num_2, inp_num_3)
			reg_mean = mean(num_vec)
			start_vec <- append(start_vec, reg_mean)
		}
		start_g <- which.max(start_vec)
		start_mat <- c(start_g, start_vec[[start_g]])
	}
	g1_start = get_ord(input)
}