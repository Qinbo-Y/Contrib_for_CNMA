####Shortest-path for CNMA
library(lpSolve)
library(MASS)

# ---- helper functions -------------------------------------------------
generate_comparisons <- function(Xa, comps, inactive) {
  comparisons <- character(nrow(Xa))
  
  for (i in seq_len(nrow(Xa))) {
    row <- Xa[i, ]
    
    left <- comps[which(row == 1)]
    right <- comps[which(row == -1)]
    
    left_str <- if (length(left) > 0) paste(left, collapse = "+") else inactive
    right_str <- if (length(right) > 0) paste(right, collapse = "+") else inactive
    
    comparisons[i] <- paste0(left_str, ":", right_str)
  }
  
  return(comparisons)
}
merge_duplicates <- function(Xa, Ha) {
  Xa_str <- apply(Xa, 1, paste, collapse = "_") 
  unique_rows <- !duplicated(Xa_str)              
  keep_idx <- which(unique_rows)
  map_to_unique <- match(Xa_str, Xa_str[keep_idx])
  
  Xa_new <- Xa[keep_idx, , drop = FALSE]
  Ha_new <- matrix(0, nrow = nrow(Ha), ncol = length(keep_idx))
  
  for (col_idx in seq_len(ncol(Ha))) {
    unique_col <- map_to_unique[col_idx]
    Ha_new[, unique_col] <- Ha_new[, unique_col] + Ha[, col_idx]
  }
  
  return(list(Xa = Xa_new, Ha = Ha_new))
}
solve_system <- function(X, e, tolerance = 1e-3) {
  
  m <- nrow(X)
  n <- ncol(X)
  Xt <- t(X)
  f.obj <- rep(0, m)
  f.con <- Xt
  f.rhs <- as.numeric(e)
  f.dir <- rep("=", n)
  
  result <- tryCatch({
    sol <- lpSolve::lp(
      direction = "min",
      objective.in = f.obj,
      const.mat = f.con,
      const.dir = f.dir,
      const.rhs = f.rhs,
      all.int = FALSE,    
      all.bin = FALSE,    
      compute.sens = FALSE
    )
    
    if (sol$status != 0) {
      return(NULL)  
    }
    
    # Coefficients for each path
    b <- sol$solution
    
    # Check whether the coefficients are non-negative
    if (any(b < -tolerance)) {
      return(NULL)
    }
    
    if (max(abs(Xt %*% b - e)) > tolerance) {
      return(NULL)
    }
    
    return(b)
  }, error = function(e) {
    return(NULL)
  })
  
  return(result)
}
component_analysis <- function(Xa, Ha, s, kmax = 4) {
  m <- nrow(Xa)
  n <- ncol(Xa)
  # Entry for each direct comparison
  edge_values <- abs(Ha[s + 1, ])
  edges <- as.list(edge_values)
  names(edges) <- as.character(0:(length(edge_values) - 1))
  stream <- list()
  
  must_have_rows <- which(Xa[, s + 1] != 0) - 1
  # identify from the shortest path
  kmin <- 1
  repeat {
    found_path <- FALSE
    
    for (k in kmin:kmax) {
      if (length(edges) < k) break
      
      combs <- list()
      utils::combn(names(edges), k, FUN = function(comb) {
        comb_int <- as.integer(comb)
        

        if (!any(comb_int %in% must_have_rows)) return(NULL)
        
        flow_sum <- sum(sapply(comb, function(i) edges[[i]]))
        combs[[length(combs) + 1]] <<- list(comb = comb, flow = flow_sum)
        return(NULL)
      }, simplify = FALSE)
      
      for (entry in combs) {
        comb_idx_chr <- entry$comb
        if (!all(comb_idx_chr %in% names(edges))) next
        
        comb_idx <- as.integer(comb_idx_chr)
        X_comb <- Xa[comb_idx + 1, , drop = FALSE]
        
        flip_mask <- Ha[s + 1, comb_idx + 1] < 0
        if (any(flip_mask)) {
          X_comb[flip_mask, ] <- -X_comb[flip_mask, ]
        }
        
        e_s <- rep(0, n)
        e_s[s + 1] <- 1
        
        sol <- tryCatch(
          solve_system(X_comb, e_s),
          error = function(e) NULL
        )
        if (is.null(sol)) next
        
        flows <- sapply(comb_idx_chr, function(i) edges[[i]])
        weighted_flows <- flows / sol
        min_flow <- min(flows)
        phi <- min(weighted_flows)
        
        for (i_chr in comb_idx_chr) {
          new_flow <- edges[[i_chr]] - min_flow
          if (new_flow < 1e-4) {
            edges[[i_chr]] <- NULL
          } else {
            edges[[i_chr]] <- new_flow
          }
        }
        
        stream[[length(stream) + 1]] <- list(sol = sol, pi = comb_idx, phi = phi)
        

        kmin <- max(kmin, k)
        found_path <- TRUE
        break  
      }
      
      if (found_path) break
    }
    
    if (!found_path) break  
  }
  
  # Identify contribution to each edge
  proportion <- numeric(length(edge_values))
  for (i in seq_along(edge_values) - 1) {
    for (entry in stream) {
      sol <- entry$sol
      pi <- entry$pi
      phi <- entry$phi
      if (i %in% pi) {
        pos_i <- which(pi == i)
        weight_sol <- sol / sum(sol)
        proportion[i + 1] <- proportion[i + 1] + phi * weight_sol[pos_i]
      }
    }
  }
  proportion <- round(proportion, 5)
  
  # Residual allocation
  leftover_edges <- list()
  if (sum(proportion) < 0.999) {
    if (length(edges) > 0) {
      leftover_sum <- sum(unlist(edges))
      remaining_weight <- 1 - sum(proportion)
      
      leftover_edges <- list()
      for (i_chr in names(edges)) {
        i <- as.integer(i_chr)
        flow <- edges[[i_chr]]
        contribution <- flow / leftover_sum * remaining_weight
        proportion[i + 1] <- proportion[i + 1] + contribution
        leftover_edges[[i_chr]] <- contribution
      }
    }
  }
  proportion <- round(proportion, 5)
  
  return(list(
    stream = stream,
    proportion = proportion,
    leftover_edges = leftover_edges
  ))
}
cnma_contrib <- function(cnma, component, inactive, random = TRUE, fixed = NULL, kmax = 4) {
  if (!is.null(fixed)) random <- !fixed
  seTE <- if (random) cnma$seTE.adj.random else cnma$seTE.adj.common
  Xa <- cnma$X.matrix
  comps <- cnma$comps
  weight <- 1 / (seTE^2)
  W <- diag(weight)
  
  # ---- Hat matrix ----
  pseudo_inverse <- MASS::ginv(t(Xa) %*% W %*% Xa)
  Ha <- pseudo_inverse %*% t(Xa) %*% W
  
  # ---- Original comparisons ----
  comparisons0 <- generate_comparisons(Xa, comps, inactive)
  
  # ---- CNMA comparisons ----
  merged <- merge_duplicates(Xa, Ha)
  Xa_m <- merged$Xa
  Ha_m <- merged$Ha
  comparisons1 <- generate_comparisons(Xa_m, comps, inactive)
  
  # ---- Component position ----
  s <- match(component, comps) - 1
  if (is.na(s)) stop("Component not found in comps vector.")
  
  # ---- Pseudo path analysis ----
  result <- component_analysis(Xa_m, Ha_m, s)
  stream <- result$stream
  proportion <- result$proportion
  
  # ---- Proportion table ----
  prop_df_filtered <- data.frame(
    Edge = comparisons1,
    Proportion = proportion,
    stringsAsFactors = FALSE
  )
  prop_df_filtered <- prop_df_filtered[prop_df_filtered$Proportion > 1e-4, ]
  
  # ---- Comparisons table ----
  treat_comparisons <- paste(cnma$treat1, cnma$treat2, sep = ":")
  study_comparisons <- data.frame(
    studyID = cnma$studlab,
    study_comparison = treat_comparisons,
    cnma_comparison = comparisons0,
    seTE_adj = seTE,
    stringsAsFactors = FALSE
  )
  
  indicator <- apply(Xa, 1, function(row) {
    col_val <- row[s + 1]
    if (!(col_val == 1 || col_val == -1)) return(FALSE)
    other_cols <- row[-(s + 1)]
    all(other_cols == 0)
  })
  study_comparisons$integrated_direct_estimate <- as.integer(indicator)
  study_comparisons <- study_comparisons[order(-study_comparisons$integrated_direct_estimate), ]
  

  stream_in_letters <- list()
  for (entry in stream) {
    sol <- entry$sol
    comb_idx <- entry$pi
    phi <- entry$phi
    

    paths <- character()
    for (i in comb_idx) {
      comp <- comparisons1[i + 1]
      parts <- strsplit(comp, ":")[[1]]
      left <- parts[1]; right <- parts[2]
      if (Ha[s + 1, i + 1] < 0) {
        path_str <- paste0('"', right, ":", left, '"')
      } else {
        path_str <- paste0('"', comp, '"')
      }
      paths <- c(paths, path_str)
    }
    

    expr <- paste(paths, collapse = " + ")
    

    b_vec <- paste(round(sol, 4), collapse = ", ")
    

    estimate_type <- if (length(comb_idx) == 1) {
      "Integrated direct estimate"
    } else {
      "Additive indirect estimate"
    }
    

    stream_in_letters[[length(stream_in_letters) + 1]] <- data.frame(
      Path = expr,
      Path_Length = length(comb_idx),
      Estimate_Type = estimate_type,
      b_vector = b_vec,
      Contrib_for_path = round(phi, 4),
      stringsAsFactors = FALSE
    )
  }
  

  Paths_df <- do.call(rbind, stream_in_letters)
  # ---- Warning if total flow is too low ----
  total_phi <- sum(Paths_df$Contrib_for_path, na.rm = TRUE)
  if (total_phi < 0.5) {
    message(sprintf(
      "⚠️ The total contribution of all pseudo paths under kmax = %d is only %.3f (< 0.5). The component '%s' may not be fully captured.",
      kmax, total_phi, component
    ))
  }
  
  # ---- Final output ----
  out <- list(
    Paths = Paths_df,
    Contributions_for_edges = prop_df_filtered,
    Comparisons_in_CNMA = study_comparisons
  )
  
  return(out)
}





# Function cnma_contrib (x, "component", "inactive", random, kmax)
#
# Arguments:
#   x: An object of 'discomb' or 'netcomb', available for both additive and interaction models..
#   component: A character string defining the target component whose contribution is to be analyzed.
#   inactive: A character string defining the inactive (reference) treatment component.
#   random：A logical indicating whether a random effects network meta-analysis should be conducted.
#   kmax: Maximum length of pseudo paths to consider. Default is 4.
#
# Returns:
#   A list with:
#   - Paths: pseudo paths and their contributions (φ) to the target component.
#   - Contributions_for_edges: A data frame of proportional contributions of each edge.
#   - Comparisons_in_CNMA: A data frame showing the original study-level comparisons, 
#     their CNMA model representation, and whether they contributed to the integrated direct estimate.



#Hypothetical data-------------------------------------------------
library(netmeta)
Study_ID <- c(1:8)
treat1 <- c("A","B","A","A+B", "A","C","A+C", "A+C")
treat2 <- c("Placebo","Placebo","B","Placebo","C","Placebo","Placebo","C")
TE <- c(-0.43, -0.39, -0.12, -0.82, 0.26, -0.51, -0.22, -0.19)
seTE <- c(0.34, 0.32, 0.33, 0.18, 0.23, 0.46, 0.29, 0.45)


#Run hypothetical data-------------------------------------------------
#NMA 
fit <- netmeta(TE, seTE, treat1, treat2, sm = "OR", 
               studlab = Study_ID, random = FALSE)
#additivel CNMA 
nc0 <- netcomb (fit, inactive = "Placebo")
# Edge combinations and proportional contributons
contrib_A <- cnma_contrib (nc0, component = "A", 
                         inactive = "Placebo", random = FALSE)


Paths_A <- contrib_A$Paths
Edges_A <- contrib_A$Contributions_for_edges
comparisons_A <- contrib_A$Comparisons_in_CNMA

