
Pretraitement_base_R <- function(bool = TRUE) {
  library(readxl)
  ######################################
  # CHARGEMENT DES DONNÉES                                   
  ######################################
  #
  Data_liege <- read_excel("data_n.xlsx")
  Data_vin <- read_excel("matrice_vin.xlsx")
  names(Data_vin)[names(Data_vin) == "Mass"] <- "mz" 
  Data_vin$mz <- as.numeric(Data_vin$mz)
  
  ############################################
  # MATCHING ET FILTRAGE DES M/Z                             
  ############################################
  #
  mz_liege_initial <- names(Data_liege)[-(1:3)]

  ppm_tol <- function(mz, ppm = 5) {
    tol <- mz * ppm / 1e6
    c(mz - tol, mz + tol)
  }
  
  # Le matching.
  match_global <- vapply(as.numeric(mz_liege_initial), function(mz_target) {
      bnds <- ppm_tol(mz_target, 5)
      potential_matches <- which(Data_vin$mz >= bnds[1] & Data_vin$mz <= bnds[2])
      if (length(potential_matches) == 0) return(NA_integer_)
      potential_matches[which.min(abs(Data_vin$mz[potential_matches] - mz_target))]
                      },integer(1))
  
  mz_avec_match <- mz_liege_initial[!is.na(match_global)]
  nbr_ech_int <- colSums(Data_liege[,mz_avec_match] > 0, na.rm = TRUE)
  int_cols_filt <- mz_avec_match[nbr_ech_int>=10]
  
  if(bool) cat("Nombre de m/z conservés après double filtrage :",length(int_cols_filt),"\n")
  ####################################################
  # PRÉPARATION DES TABLES LIÈGE ET VIN                      
  ####################################################
  #
  cols_a_garder <- c("Millesime","Masse","Position",int_cols_filt)
  Data_filt <- Data_liege[,cols_a_garder]
  Data_filt[,int_cols_filt] <- sweep(Data_filt[,int_cols_filt],1,Data_filt$Masse, FUN = "/")
  Data_filt[,int_cols_filt] <- log(Data_filt[,int_cols_filt] + 1)
  
  match <- match_global[mz_liege_initial %in% int_cols_filt]
  vin_T <- t(Data_vin[match, setdiff(colnames(Data_vin), "mz")])
  colnames(vin_T) <- int_cols_filt
  
  vin_r <- data.frame(Millesime = rownames(vin_T),
                          Masse = NA_real_,
                          Position = 6,
                          log(vin_T+1), # Transformation log appliquée directement
                          check.names = FALSE, # Empêche R de modifier les noms de colonnes (ex: avec des "X")
                          row.names = NULL     # Supprime les anciens noms de lignes
  )
  
  data_final <- rbind(Data_filt, vin_r)
  
  return(data_final[, cols_a_garder]) 
}





FonctionPlsDa <- function(bool=TRUE){
  #PACKAGE ##########
  library(mixOmics) #
  #DONNEES ############################
  data_filt <- Pretraitement_base_R() #
  X <- data_filt[,-c(1:3)]            #
  Y <- as.factor(data_filt$Position)  #
  #ANALYSE ###############################################################
  plsda_model <- mixOmics::plsda(X,Y,ncomp=2,scale = TRUE)               #
  vip_scores <- vip(plsda_model)                                         #
  vip_df <- data.frame(mz = rownames(vip_scores),VIP = vip_scores[, 1])  #
  vip_df_sorted <- vip_df[order(vip_df$VIP, decreasing = TRUE),]         #
  
  if(bool==TRUE){
    mixOmics::plotIndiv(plsda_model, comp = c(1, 2), group = Y,title = "PLS-DA (Position)", legend = TRUE, ellipse = TRUE)
    plotVar(plsda_model, comp = c(1, 2), title = "Variables sélectionnées (comp1)")
    top_vip_n <- vip_df_sorted[1:20, ]
    print(top_vip_n)
  }
  
  return (vip_df_sorted$mz[1:20])
}


FonctionPlsDa(TRUE)
