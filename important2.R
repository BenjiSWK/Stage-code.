library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(mixOmics)


#  mds  -  tsne   -   15 16 17

#############################
# PRETRAITEMENT DES DONNEES #
#############################

Pretraitement <- function(bool=TRUE){
data <- read_excel("data_n .xlsx")
data <- data[,order(colnames(data))]
#FILTRE VIN
mz_vin <- as.numeric(read_excel("matrice_vin.xlsx")$Mass)
ppm_tol <- function(mz,ppm){
  tol <- mz*ppm/1e6
  cbind(lower = mz-tol,upper = mz+tol)
}


bounds <- ppm_tol(mz_vin, ppm = 5)
mz_data <- as.numeric(colnames(data)[-1:-3])
keep <- sapply(mz_data, function(x) any(x >= bounds[,1] & x <= bounds[,2]))
mz_migres <- mz_data[keep]

mat_liege_migres <- as.matrix(data[,(4:ncol(data))])[,keep]
data <- cbind(data[,1:3],as.data.frame(mat_liege_migres))


cat("m/z initiaux dans data_n :",length(mz_data),"\n")
cat("m/z après matching ppm :",length(mz_migres),"\n")




data <- data[data$Position %in% c(1,2,3,4,5),]
# Filtrage des m/z peu informatifs
colonne_intensité <- names(data)[-(1:3)]  
nbr_ech_int <- colSums(data[, colonne_intensité] > 0, na.rm = TRUE) 
seuil <- 10
int_cols_filt <- names(nbr_ech_int)[nbr_ech_int>=seuil]

data_filt <- data[, c("echantillon", "gramme", "Position", int_cols_filt)] 
if(bool==TRUE){
#cat("Nombre de colonnes supprimées aprés filtage cork :", length(cols_to_drop), "\n")
cat("Nombre de variables m/z avant filtrage : 12958 ", "\n")
cat("Nombre de variables m/z après filtrage  :", length(int_cols_filt), "\n")
}
# Normalisation par le poids 
data_filt[, int_cols_filt] <- sweep(data_filt[,int_cols_filt], 1, data_filt$gramme, FUN = "/")
# Transformation log des intensités
data_filt[, int_cols_filt] <- log(data_filt[, int_cols_filt] + 1)

#data_filt$Millesime <- as.numeric(sub("\\(.*", "", data_filt$echantillon))
#data_filt <- data_filt[data_filt$Millesime >= 1989, ]




# 
# 
# 
# 
# mz_cols <- names(data_filt)[-(1:3)]    # toutes les colonnes après echantillon, gramme, Position
# mz_vals <- as.numeric(mz_cols)         # convertit ces noms en vecteur numérique
# 
# # 1. Lire la matrice “vin” et renommer Mass → mz
# vin_raw <- readxl::read_excel("matrice_vin.xlsx") %>%
#   dplyr::rename(mz = Mass)
# 
# # 2. Forcer mz en numérique
# vin_raw$mz <- as.numeric(vin_raw$mz)
# 
# # 3. Déterminer les colonnes correspondant aux échantillons
# #    (tout ce qui n’est pas la colonne mz)
# sample_cols <- setdiff(colnames(vin_raw), "mz")
# 
# # 4. Fonction pour bornes ppm
# ppm_tol <- function(mz, ppm = 5){
#   tol <- mz * ppm/1e6
#   c(mz - tol, mz + tol)
# }
# 
# # 5. Pour chaque mz_vals, trouver l’indice dans vin_raw$mz
# match_idx <- vapply(mz_vals, function(mz) {
#   bnds <- ppm_tol(mz, 5)
#   idxs <- which(vin_raw$mz >= bnds[1] & vin_raw$mz <= bnds[2])
#   if(length(idxs) == 0) return(NA_integer_)
#   # si plusieurs, prendre le plus proche
#   idxs[which.min(abs(vin_raw$mz[idxs] - mz))]
# }, integer(1))
# 
# # 6. Extraire la sous‐matrice (lignes = mz_vals, colonnes = samples)
# matched_mat <- vin_raw[match_idx, sample_cols, drop = FALSE]
# 
# # 7. Donner à chaque ligne le nom de mz_vals
# rownames(matched_mat) <- as.character(mz_vals)
# 
# # 8. Transposer pour que lignes = échantillon, colonnes = mz
# vin_wide <- as.data.frame(t(matched_mat), stringsAsFactors = FALSE)
# vin_wide$echantillon <- rownames(vin_wide)
# 
# # 9. Ajouter gramme = NA et Position = 6
# vin_wide$gramme   <- NA_real_
# vin_wide$Position <- 6
# 
# # 9. Après avoir construit vin_wide et ajouté Position=6 :
# vin_wide[, int_cols_filt] <- log(vin_wide[, int_cols_filt] + 1)
# 
# # 10. Puis on réordonne et on fusionne
# vin_wide <- vin_wide[, c("echantillon","gramme","Position", int_cols_filt)]
# data_final <- rbind(data_filt, vin_wide)

#data_filt$Millesime <- as.integer(
#  sub("\\(.*", "", data_filt$echantillon)
#)


return (data_filt)
}














                         ###################
              ############ Fonction PLS-DA ############
                         ###################



FonctionPlsDa <- function(bool=TRUE){
  library(mixOmics)
  library(pls)
  data_filt <- Pretraitement(FALSE)
  
  # Préparation des données
  X <- data_filt[,-c(1:3)]
  #Y <- factor(ifelse(data_filt$Position == 1,"Pos1",ifelse(data_filt$Position == 2, "Pos2", "Pos345")),levels = c("Pos1","Pos2","Pos345"))
  Y <- as.factor(data_filt$Position)
  
  # PLS-DA
  plsda_model <- mixOmics::plsda(X,Y,ncomp=2)
  
  # VIP scores
  score_vip <- vip(plsda_model)
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1), ]
  
  if(bool==TRUE){
  # Affichage des individus
  mixOmics::plotIndiv(plsda_model, comp = c(1, 2), group = Y,title = "PLS-DA (Position)", legend = TRUE, ellipse = TRUE)
  plotVar(plsda_model, comp = c(1, 2), title = "Variables sélectionnées (comp1)")
  top_vip_n <- vip_df_sorted[1:20, ]
  print(top_vip_n)
  }
  return (vip_df_sorted$mz[1:20])
}



                   ######################
      ############## FONCTION AFFICHAGE #############
                   ######################

AffichageBoiteMoustache <- function(top_mz){
  
  data_filt<-Pretraitement(FALSE)
  
  # Format long pour ggplot
  library(reshape2)
  data_plot <- melt(data_filt[, c("Position", top_mz)],
                    id.vars = "Position",
                    variable.name = "m_z",
                    value.name = "intensité")
  
  # Tracer les profils selon la Position
  ggplot(data_plot, aes(x = as.factor(Position), y = intensité, fill = Position)) +
    geom_boxplot() + facet_wrap(~ m_z, scales = "free_y") + theme_minimal() +
    labs(title = "Profils des 20 m/z les plus discriminants (VIP)",x = "Position", y = "Intensité")
  
}




                        #################
           ############## PROFILS MOYEN #############
                        #################


AffichageProfilsMoyen <- function(top_mz){
  library(dplyr)
  data_filt<-Pretraitement(FALSE)
  
  melted_data <- reshape2::melt(data_filt[,c("Position",top_mz)],id.vars = "Position",variable.name = "m_z",value.name = "intensité")
  
  # Moyennes et écart-type
  summary_df <- melted_data %>% group_by(m_z, Position) %>% summarise(mean_intensity = mean(intensité, na.rm = TRUE),
                sd = sd(intensité, na.rm = TRUE), n = n(), se = sd / sqrt(n), .groups = "drop")
  
  # Tracé des profils moyens avec SE
  ggplot(summary_df, aes(x = as.numeric(Position),y = mean_intensity,group = m_z,color = m_z)) +
    geom_line() + geom_point() + geom_errorbar(aes(ymin = mean_intensity - se,ymax = mean_intensity + se),width = 0.1) +
    facet_wrap(~ m_z, scales = "free_y") + theme_minimal() +labs(title = "Profil moyen d’intensité selon la profondeur",
    x = "Position dans le bouchon (1 = contact vin)", y = "Intensité moyenne (normalisée)")
  
}


                            ####################
               ############## FONCTION HEATMAP #############
                            ####################



AffichageHeatmap <- function(top_mz){
  data_filt<-Pretraitement(FALSE)
  
  # Réordonner les échantillons selon Position
  ord <- order(data_filt$Position)
  X_heat_ord <- data_filt[ord, top_mz]
  rownames(X_heat_ord) <- paste(data_filt$echantillon[ord], data_filt$Position[ord], sep = "_")
  
  # Créer l'annotation colonne
  annotation_col <- data.frame(Position = as.factor(data_filt$Position[ord]))
  rownames(annotation_col) <- rownames(X_heat_ord)
  
  # Transposer pour que les m/z soient en lignes
  pheatmap(t(X_heat_ord),annotation_col = annotation_col,scale = "row",cluster_rows = TRUE,cluster_cols = FALSE,show_colnames = FALSE)
}


                           #####################
              ############## FONCTION SPEARMAN #############
                           #####################


SpearmanMigrateurs <- function() {
  data_filt <-  Pretraitement(FALSE)
  seuil <- -0.5
  mz_names <- colnames(data_filt)[-(1:3)]
  res <- data.frame(mz = mz_names, cor = NA)
  
  for (i in seq_along(mz_names)) {
    intensite <- data_filt[[mz_names[i]]]
    pos <- as.numeric(data_filt$Position)
    res$cor[i] <- suppressWarnings(cor(intensite, pos, method = "spearman", use = "complete.obs"))
  }
  
  mz_migrateurs <- res$mz[res$cor < seuil]
  print(res[res$cor < seuil, ])  
  return(mz_migrateurs)
}



                           ##########################
              ############## FONCTION TRANSPOSITION #############
                           ##########################




FonctionTranspositionMatrice <- function(){
  #PACKAGES###########
  library(FactoMineR)#
  library(factoextra)#
  library(mixOmics)  #
  #DONNEES#############################
  data <- Pretraitement(FALSE)        #
  mat <- data.matrix(data[,-c(1:3)]) # 
  X <- t(mat)                         #
  Y <- as.factor(t(data$echantillon)) #
  #ANALYSE#################################################################################################################################################
  res_pca <- PCA(X,scale.unit = TRUE,graph = FALSE)                                                                                                       
  print(fviz_pca_ind(res_pca, label = "none",mean.point = FALSE)+labs(title ="Projection des échantillons"))                                                    
  print(fviz_pca_var(res_pca,axes = c(1, 2),habillage = Y,label="none",geom = c("point","text")) +labs(title = "PCA (transposée) : projection des échantillons"))
}




                       ###################
          ############## PLSDA_Millesime #############
                       ###################



PLSDA_Millesime <- function(){
  #PACKAGES#########
  library(mixOmics)#
  #DONNEES##################################
  data_filt  <- Pretraitement(FALSE)       #
  X <- as.matrix(data_filt[ , -c(1:3)])    #
  Y_ech <- as.factor(data_filt$echantillon)#
  Y_pos <- as.factor(data_filt$Position)   #  
  #ANALYSE#########################################################################################################################
  plsda_model  <- plsda(X, Y_ech, ncomp = 2)                                                                                     #
  plotIndiv(plsda_model,comp = c(1, 2),group = Y_ech,legend = TRUE,title = "PLS-DA : point labels = Position")#
  plotVar(plsda_model, comp = c(1, 2), title = "Variables sélectionnées (comp1)")
  score_vip <- vip(plsda_model)                                                                                                  #
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])                                                    #
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1), ]                                                                    #
  return(vip_df_sorted$mz[1:20])                                                                                                  #
}

                      ####################
         ############## sPLSDA_Millesime #############
                      ####################


SPLSDA_Millesime <- function(bool=TRUE){
  #PACKAGE############
  library(mixOmics)  #
  library(factoextra)#
  #DONNEES#######################################
  data_filt  <- Pretraitement(FALSE)            #
  X          <- as.matrix(data_filt[ , -c(1:3)])#
  Y_ech      <- as.factor(data_filt$echantillon)#
  Y_pos      <- as.factor(data_filt$Position)   #
  #ANALYSE#########################################################################################################################
  splsda_model  <- splsda(X, Y_ech, ncomp = 2,keepX = c(40, 40))
  if(bool==TRUE){
  mixOmics::plotIndiv(splsda_model,comp = c(1, 2),group = Y_ech,ind.names = Y_pos,legend = TRUE,title = "sPLS-DA : point labels = Position")#
  print(plot_vars(splsda_model))
  }
  score_vip <- vip(splsda_model)                                                                                                  #
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])                                                    #
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1), ]                                                                    #
  return(vip_df_sorted$mz[1:20])                                                                                                  #
}

                       ############
          ############## Plot var #############
                       ############


plot_vars <- function(model) {
  library(ggplot2)
  library(ggrepel)
  
  # Récupération des coordonnées
  resVar <- plotVar(model, comp = c(1,2), plot = FALSE)
  df <- data.frame(comp1 = resVar$x,comp2 = resVar$y,var   = resVar$names)
  ggplot(df, aes(x = comp1, y = comp2)) +
  geom_path(data = data.frame(theta = seq(0, 2*pi, length.out = 200)),aes(x = cos(theta), y = sin(theta)),inherit.aes = FALSE,linetype = "dashed",colour = "grey80") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 1) +
  geom_text_repel(aes(label = var),size = 2,max.overlaps = Inf,box.padding = 0.5,point.padding = 0.3,segment.size  = 0.2,segment.color = "grey50") +
  coord_equal() + labs(title = "Variables sélectionnées (comp1 vs comp2)",x = paste("Composante", c(1)),y = paste("Composante", c(2))) +
  theme_minimal()
}




                            ####################
               ############## FONCTION sPLS-DA #############
                            ####################




FonctionSplsDa <- function(bool = TRUE) {
  library(mixOmics)
  data_filt <- Pretraitement(FALSE)
  X <- data_filt[, -c(1:3)]
  #Y <- factor(ifelse(data_filt$Position == 1,"Pos1",ifelse(data_filt$Position == 2, "Pos2", "Pos345")),levels = c("Pos1","Pos2","Pos345"))
  Y <- as.factor(data_filt$Position)
  
  keepX = c(20, 20)

  splsda_model <- mixOmics::splsda(X, Y, ncomp = 2, keepX = c(20, 20))
  score_vip <- vip(splsda_model)
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1),]
  
  if (bool) {
    mixOmics::plotIndiv(splsda_model,comp = c(1,2),group = Y,
    title = "sPLS-DA (Position)",legend = TRUE,ellipse = TRUE)
    print(plot_vars(splsda_model))
    top_vip_n <- vip_df_sorted[1:20, ]
    print(top_vip_n)
  }

  return(vip_df_sorted$mz[1:20])
}

                 ###########
          ######## PCA VIN #######
                 ###########

PCA_VIN <- function(){
  library(FactoMineR)  
  library(factoextra)
  data <- readxl::read_excel("matrice_vin.xlsx")
  X <- t(data[,-c(1)])
  Y <- as.factor(colnames(data[,-c(1)]))
  res_pca <- PCA(X,scale.unit = FALSE,graph = FALSE)
  print(fviz_pca_ind(res_pca,habillage = Y, label = "none",mean.point = FALSE,addEllipses = FALSE) +labs(title ="Projection des échantillons"))
  contrib_total <- rowSums(get_pca_var(res_pca)$contrib[,1:2])
  top10_total <- names(sort(contrib_total, decreasing = TRUE))[1:10]
  fviz_pca_var(res_pca, select.var = list(name = top10_total),col.var = "contrib")+labs(title ="TOP 10 contribution des variables")
}

            #######
     ######## PCA #######
            #######

PCA <- function(){
  library(FactoMineR)  
  library(factoextra)
  data <- Pretraitement(FALSE)
  X <- data[, -c(1:3)]
  #Y <- factor(ifelse(data_filt$Position == 1,"Pos1",ifelse(data_filt$Position == 2, "Pos2", "Pos345")),levels = c("Pos1","Pos2","Pos345"))
  Y <- as.factor(data$Position)
  Y_ech <- as.factor(data$Millesime)
  res_pca <- PCA(X,scale.unit = TRUE,graph = FALSE)
  print(fviz_pca_ind(res_pca,habillage = Y, label = "none",mean.point = FALSE,addEllipses = TRUE,pointshape  = 19) +labs(title ="Projection des échantillons"))
  print(fviz_pca_ind(res_pca,habillage = Y_ech, label = "none",mean.point = FALSE,addEllipses =TRUE,pointshape  = 19) +labs(title ="Projection des échantillons"))
  
  contrib_total <- rowSums(get_pca_var(res_pca)$contrib[,1:2])
  top10_total <- names(sort(contrib_total, decreasing = TRUE))[1:10]
  fviz_pca_var(res_pca, select.var = list(name = top10_total),col.var = "contrib")+labs(title ="TOP 10 contribution des variables")
}

             ###################
      ######## Foret Aleatoire #######
             ###################

ForetAleatoire <- function(){
  # PACKAGES ############
  library(randomForest) #
  library(dplyr)        #
  
  df <- Pretraitement(FALSE)
  Y <- as.factor(df$Position)
  X <- df[,-c(1:3)]
  
  rf_mod <- randomForest(x=X,y=Y,ntree=500,importance=TRUE)
  
  cat("Erreur OOB (overall) :", round(rf_mod$err.rate[500, "OOB"], 3), "\n")
  
  imp <- importance(rf_mod, type = 2)  # type=2 -> Gini
  imp_df <- data.frame(mz=rownames(imp),MeanDecreaseGini=imp[,"MeanDecreaseGini"]) 
  imp_df <- imp_df[order(imp_df$MeanDecreaseGini, decreasing = TRUE),]
  top10_rf <- head(imp_df, 10)
  print(top10_rf)
  
  varImpPlot(rf_mod, n.var = 10, main = "Top 10 m/z - Random Forest (Gini)")
  
  err_mat <- rf_mod$err.rate
  print(err_mat)
  
  conf_mat <- rf_mod$confusion
  print(conf_mat)
}

             #######
      ######## HCA #######
             #######


HCA <- function(){
  library(pheatmap)
  
  df <- Pretraitement(FALSE)
  X <- as.matrix(df[, -(1:3)])     
  ids <- paste0(df$echantillon, "_pos", df$Position)
  rownames(X) <- ids
  
  X[!is.finite(X)] <- NA
  X <- X[complete.cases(X), , drop = FALSE]
  vars <- apply(X, 2, var, na.rm = TRUE)
  X <- X[, vars > 0, drop = FALSE]
  
  dist_mat <- dist(X, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  plot(hc,labels = rownames(X),main="Dendrogramme HCA des échantillons",xlab= "",sub = "",cex = 0.7)
  
  ord <- order(df$Position)
  X_ord <- X[ord, , drop = FALSE]
  ids_ord <- ids[ord]
  
  annotation_col <- data.frame(Position = factor(df$Position[ord]))
  rownames(annotation_col) <- ids_ord
  
  pheatmap(t(X_ord),annotation_col = annotation_col,scale = "row",cluster_rows = TRUE, cluster_cols     = FALSE,   # ne pas recluster les colonnes
           show_rownames = FALSE,show_colnames = TRUE,main = "Heatmap HCA — Position triée")
}






  ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    
    ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ###






                    ###########
          ########### PLS-DA  ###########
                    ###########

FonctionPlsDa(TRUE)
AffichageBoiteMoustache(FonctionPlsDa(FALSE))
AffichageProfilsMoyen(FonctionPlsDa(FALSE))
AffichageHeatmap(FonctionPlsDa(FALSE))


                    ###########
         ############ sPLS-DA ############
                    ###########

FonctionSplsDa(TRUE)
AffichageBoiteMoustache(FonctionSplsDa(FALSE))
AffichageProfilsMoyen(FonctionSplsDa(FALSE))
AffichageHeatmap(FonctionSplsDa(FALSE))


                    #########################
         ############ TRANSPOSITION MATRICE ############
                    #########################

FonctionTranspositionMatrice()

                   ########################
         ########### PLS-DA PAR MILLESIME ###########
                   ########################

PLSDA_Millesime()
AffichageBoiteMoustache(PLSDA_Millesime())
AffichageProfilsMoyen(PLSDA_Millesime())
AffichageHeatmap(PLSDA_Millesime())

                   #########################
         ########### sPLS-DA PAR MILLESIME ###########
                   #########################

SPLSDA_Millesime(TRUE)
AffichageBoiteMoustache(SPLSDA_Millesime(FALSE))
AffichageProfilsMoyen(SPLSDA_Millesime(FALSE))
AffichageHeatmap(SPLSDA_Millesime(FALSE))

                     #######
              ######## PCA #######
                     #######

PCA_VIN()

               ###################
        ######## Foret Aleatoire #######
               ###################

ForetAleatoire()




















# 1. Chargement des packages
library(ggplot2)

# 2. Chargement de vos données prétraitées
#    On part du résultat de Pretraitement(FALSE) ou data_final s'il inclut la position 6
df <- Pretraitement(FALSE)      # ou data_filt si vous ne voulez pas la position 6

# 3. Identification des colonnes d'intensité m/z
mz_cols <- names(df)[!(names(df) %in% c("echantillon","gramme","Position"))]

cor_mat  <- cor(t(df[, mz_cols]), use = "pairwise.complete.obs")
dist_mat <- as.dist(1 - cor_mat)

# 5. Exécution de la MDS (classique)
mds_res <- cmdscale(dist_mat, k = 2, eig = TRUE)

# 6. Construction d'un data.frame pour le plot
mds_df <- data.frame(
  MDS1      = mds_res$points[, 1],
  MDS2      = mds_res$points[, 2],
  Position  = factor(df$Position)
)

# 7. Visualisation avec ggplot2
ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Position)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title    = "MDS des profils m/z",
    x        = "Dimension 1",
    y        = "Dimension 2",
    color    = "Position",
    shape    = "Millésime"
  ) +
  theme(
    legend.position = "right"
  )

