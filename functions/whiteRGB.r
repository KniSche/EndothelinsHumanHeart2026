  whiteRGB <<- function(r_values, g_values, b_values) {
  #    """
  #    Alternative implementation: Each channel represents intensity from white
  #    Higher values move that channel away from white (1) towards 0
  #    """
      
      # Validate inputs (same as main function)
      if (!is.numeric(r_values) || !is.numeric(g_values) || !is.numeric(b_values)) {
        stop("All inputs must be numeric vectors")
      }
      if (length(r_values) != length(g_values) || length(r_values) != length(b_values)) {
        stop("All input vectors must have the same length")
      }
      if (any(r_values < 0 | r_values > 1) || 
          any(g_values < 0 | g_values > 1) || 
          any(b_values < 0 | b_values > 1)) {
        stop("All values must be between 0 and 1")
      }
      
      # Simple approach: intensity directly controls how much to reduce from white
      final_r <- 1 - (1 - r_values) * (g_values + b_values > 0)
      final_g <- 1 - (1 - g_values) * (r_values + b_values > 0)  
      final_b <- 1 - (1 - b_values) * (r_values + g_values > 0)
      
      # Alternative: more intuitive approach
      # Each channel intensity reduces the OTHER channels
      final_r <- pmax(r_values, 1 - g_values - b_values)
      final_g <- pmax(g_values, 1 - r_values - b_values)
      final_b <- pmax(b_values, 1 - r_values - g_values)
      
      # Ensure values stay within 0-1 range
      final_r <- pmax(0, pmin(1, final_r))
      final_g <- pmax(0, pmin(1, final_g))
      final_b <- pmax(0, pmin(1, final_b))
      
      colors <- rgb(final_r, final_g, final_b)
      return(colors)
    }   
     
