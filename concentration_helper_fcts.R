#============== Concentration measures ===================#
library("MASS")
library("stats")


#======== Custom functions ============#

discr_distr_fun <- function(vec){
  order_vec <- sort(vec, decreasing = F)
  count.df <- aggregate(order_vec, by = list(order_vec), FUN = length)
  colnames(count.df) <- c("Y", "Freq")
  count.df$Prob <- count.df$Freq/sum(count.df$Freq)
  return(count.df)
}

cum_Exp_fun <- function(mat){
  out <- cumsum(x = mat[,1]*mat[,2])
  return(out)
}

Gini_coeff_fun <- function(vec){
  if(sum(abs(vec))==0){
    G <- 0
  }else{
    distr_df <- discr_distr_fun(vec)
    mat <- distr_df[, c("Y", "Prob")]  
    S_i <- cum_Exp_fun(mat)
    sum_shifted_S_i <- c(0,S_i) + c(S_i,0) 
    sum_shifted_S_i <- sum_shifted_S_i[-length(sum_shifted_S_i)]
    G <- 1-(sum(mat$Prob*sum_shifted_S_i))/S_i[length(S_i)]  
  }
  return(G)
}



# Herfindhal-Hirschman Index
HH_index_fun <- function(x){
  cum <- sum(x)
  # Normalize vector
  if(cum == 0){
    stop("No concentration: Zero vector as input.")
  }else{
    norm_x <- x/sum(x)
  }
  return(sum(norm_x^2))
}


### Lorenz Curve
discr_Lorenz_curve_fun <- function(probas, y, plot_curve = TRUE){
  # input checks
  if(length(probas) != length(y)){
    stop("Input vectors are not of the same length.")
    }
  if(sum((y[-1]-y[-length(y)])<0)>0){
    stop("Vector y is not strictly increasing.")
  }
  
  # Cum fun
  cumDistr <- c(0, cumsum(probas))
  # 
  S <- cumsum(y)
  L <- c(0, S/sum(y))
  
  lorenz_curve_pts <- matrix(data = c(cumDistr, L), byrow = FALSE, ncol = 2)
  
  if(plot_curve == TRUE){
    plot(x = lorenz_curve_pts[,1], y = lorenz_curve_pts[,2], type = "l", col = "red", 
         main="Lorenz and Cumulative Sum Curves",
         ylab=""
         )
    lines(x = c(0,1), y = c(0,1))
    lines(x = seq(from = 0, to = 1, by = 1/length(y)), y = c(0, cumsum(sort(x = y, decreasing = FALSE))/sum(y)), col = "blue")
    legend("topleft",
           c("Lorenz","Cum. Sum"),
           fill=c("red","blue")
    )
     }
  
  return(lorenz_curve_pts)
}
#============== PLOTS AND OTHER ==============#

## Function to perform multiple plots
multiple_plot_fun <- function(mat){
  if(ncol(mat)>5){
    colors_vec <- sample(colours(), ncol(mat))
  }else{
    colors_vec <- c("red", "orange", "blue", "gray")    
  }
  plot(mat[,1], type = "l")
  for(i in 2:ncol(mat)){
    par(new = T)
    lines(mat[,i], col= colors_vec[(i-1)])
  }
}


#=========== Inverse quantile function ===============#
inv_quantile_fun <- function(distr, qval){
  out <- c()
  tmp_bound_min <- 0
  ordered_qval <- sort(qval, decreasing = FALSE)
  for(i in 1:length(ordered_qval)){
    tmp_qval <- ordered_qval[i]
    quant_fun <-function(u) {as.numeric(quantile(x = distr, probs = u))-tmp_qval}
    tmp_root <- uniroot(f = quant_fun,
                        lower = tmp_bound_min, upper = 1, tol = 1e-10)
    out[i] <- tmp_root$root
    tmp_bound_min <- tmp_root$root-max(c(0,tmp_root$estim.prec), na.rm = TRUE)
  }
  out <- matrix(data = c(ordered_qval, out), byrow = TRUE, nrow = 2)
  rownames(out) <- c("Prob", "Quant_Value")
  return(out)
}
