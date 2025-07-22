

Pretraitement <- function(bool=TRUE){

# PACKAGE ############################################################
library(readxl)
######################################################################

# DONNEES ############################################################
data <- read_excel("data_n.xlsx")
######################################################################

# FILTRE VIN #########################################################
# mz_vin <- as.numeric(read_excel("matrice_vin.xlsx")$Mass)
# 
# ppm_tol <- function(mz,ppm){
#   tol <- mz*ppm/1e6
#   cbind(lower = mz-tol,upper = mz+tol)
# }
# 
# bornes <- ppm_tol(mz_vin, ppm = 5)
# mz_data <- as.numeric(colnames(data)[-1:-3])
# keep <- sapply(mz_data, function(x) any(x >= bornes[,1] & x <= bornes[,2]))
# mz_migres <- mz_data[keep]
# mat_liege_migres <- as.matrix(data[,(4:ncol(data))])[,keep]
# data <- cbind(data[,1:3],as.data.frame(mat_liege_migres))
######################################################################

# FILTRAGE DES M/Z PEU INFORMATIFS ###################################
colonne_intensité <- names(data)[-(1:3)]  
nbr_ech_int <- colSums(data[, colonne_intensité] > 0, na.rm = TRUE) 
seuil <- 10
int_cols_filt <- names(nbr_ech_int)[nbr_ech_int>=seuil]
data_filt <- data[,c("Millesime", "Masse", "Position", int_cols_filt)] 
######################################################################

if(bool==TRUE){
# cat("m/z initiaux dans data_n :",length(mz_data),"\n")
# cat("m/z après matching ppm vin:",length(mz_migres),"\n")
cat("Nombre de variables m/z après filtrage des m/z peu informatifs :",length(int_cols_filt), "\n")
}

# NORMALISATION PAR LE POIDS #########################################
data_filt[, int_cols_filt] <- sweep(data_filt[,int_cols_filt], 1, data_filt$Masse, FUN = "/")
######################################################################

# LOG-TRANSFORMATION DES INTENSITES ##################################
data_filt[, int_cols_filt] <- log(data_filt[, int_cols_filt] + 1)
######################################################################

return (data_filt)
}



                         ###################
              ############ Fonction PLS-DA ############
                         ###################



FonctionPlsDa <- function(bool=TRUE){
  #PACKAGE ##########
  library(mixOmics) #
  #DONNEES ############################
  data_filt <- Pretraitement(TRUE)    #
  X <- data_filt[,-c(1:3)]            #
  Y <- as.factor(data_filt$Position)  #
  #ANALYSE ##############################################################
  plsda_model <- mixOmics::plsda(X,Y,ncomp=2,scale = TRUE)              #
  vip_scores <- vip(plsda_model)                                        #
  vip_df <- data.frame(mz = rownames(vip_scores),VIP = vip_scores[,1])  #
  vip_df_sorted <- vip_df[order(vip_df$VIP, decreasing = TRUE),]        #

  if(bool==TRUE){
    # Affichage des individus
    mixOmics::plotIndiv(plsda_model,comp = c(1,2), group = Y,title = "PLS-DA (Position)", legend = TRUE, ellipse = TRUE)
    plotVar(plsda_model, comp = c(1,2), title = "Variables sélectionnées (comp1)")
    top_vip_n <- vip_df_sorted[1:20,]
    print(top_vip_n)
  }
  return (vip_df_sorted$mz[1:20])
}



                   ##########################
      ############## Fonction Affichage BAM #############
                   ##########################



AffichageBoiteMoustache <- function(top_mz){
  #PACKAGES #########
  library(reshape2) #
  library(ggplot2)  #
  #DONNEES ########################
  data_filt<-Pretraitement(FALSE) #
  #ANALYSE ##################################################################################################################
  data_plot <- melt(data_filt[, c("Position", top_mz)],id.vars = "Position",variable.name = "m_z",value.name = "intensité") #
  #AFFICHAGE #####################################################################################################################
  ggplot(data_plot, aes(x = as.factor(Position), y = intensité, fill = Position)) +
    geom_boxplot() + facet_wrap(~ m_z, scales = "free_y") + theme_minimal() + 
    geom_jitter(width = 0.2, alpha = 0.7, color = "black") +
    labs(title = "Profils des 20 m/z les plus discriminants (VIP)",x = "Position", y = "Intensité")
  
}



AffichageBoiteMoustacheMillesime <- function(top_mz){
  #PACKAGES #########
  library(reshape2) #
  library(ggplot2)  #
  #DONNEES ####################################################################################################################
  data_filt<-Pretraitement(FALSE)                                                                                             #
  data_filt$Millesime <- as.factor(sub("\\(.*", "", data_filt$Millesime))                                                     #
  data_plot <- melt(data_filt[, c("Millesime", top_mz)],id.vars = "Millesime",variable.name = "m_z",value.name = "intensité") #
  #AFFICHAGE ##################################################################################################################
  ggplot(data_plot, aes(x = as.factor(Millesime), y = intensité, fill = Millesime)) +
    geom_boxplot() + facet_wrap(~ m_z, scales = "free_y") + theme_minimal() +
    labs(title = "Profils des 20 m/z les plus discriminants (VIP)",x = "Millesime", y = "Intensité")
}



                        ##########################
           ############## Fonction Profils Moyen #############
                        ##########################


AffichageProfilsMoyen <- function(top_mz){
  #PACKAGE ##########
  library(reshape2) #
  library(ggplot2)  #
  library(dplyr)    #
  #DONNEES ############################################################################################################################
  data_filt<-Pretraitement(FALSE)                                                            
  melted_data <- reshape2::melt(data_filt[,c("Position",top_mz)],id.vars = "Position",variable.name = "m_z",value.name = "intensité")
  summary_df <- melted_data %>% group_by(m_z,Position) %>% summarise(mean_intensity = mean(intensité,na.rm = TRUE),
                sd = sd(intensité, na.rm = TRUE), n = n(), se = sd / sqrt(n), .groups = "drop")
  #AFFICHAGE ##############################################################################################################################
  ggplot(summary_df, aes(x = as.numeric(Position),y = mean_intensity,group = m_z,color = m_z)) +
    geom_line() + geom_point() + geom_errorbar(aes(ymin = mean_intensity - se,ymax = mean_intensity + se),width = 0.1) +
    facet_wrap(~ m_z, scales = "free_y") + theme_minimal() +labs(title = "Profil moyen d’intensité selon la profondeur",
    x = "Position dans le bouchon (1 = contact vin)", y = "Intensité moyenne (normalisée)")
  
}


                            ####################
               ############## Fonction HeatMap #############
                            ####################



AffichageHeatmap <- function(top_mz){
  #PACKAGE ##########
  library(pheatmap) #
  #DONNEES ####################################################################################
  data_filt<-Pretraitement(FALSE)                                                             #
  ord <- order(data_filt$Position)                                                            #
  X_heat_ord <- data_filt[ord,top_mz]                                                         #
  rownames(X_heat_ord) <- paste(data_filt$Millesime[ord],data_filt$Position[ord],sep = "_")   #
  annotation_col <- data.frame(Position = as.factor(data_filt$Position[ord]))                 #
  rownames(annotation_col) <- rownames(X_heat_ord)                                            #
  #AFFICHAGE ##################################################################################################
  pheatmap(t(X_heat_ord),annotation_col = annotation_col,scale = "row",cluster_rows = TRUE,cluster_cols = FALSE,show_colnames = FALSE)
}



                           ##########################
              ############## Fonction Transposition #############
                           ##########################




FonctionTranspositionMatrice <- function(){
  #PACKAGES###########
  library(FactoMineR)#
  library(factoextra)#
  library(mixOmics)  #
  #DONNEES#############################
  data <- Pretraitement(FALSE)        #
  mat <- data.matrix(data[,-c(1:3)])  # 
  X <- t(mat)                         #
  Y <- as.factor(t(data$Millesime))   #
  #ANALYSE#######################################################################################################################################################
  res_pca <- FactoMineR::PCA(X,scale.unit = TRUE,graph = FALSE)                                                                                                 #            
  print(fviz_pca_ind(res_pca, label = "none",mean.point = FALSE)+labs(title ="Projection des échantillons"))                                                    #
  print(fviz_pca_var(res_pca,axes = c(1,2),habillage = Y,label="none",geom = c("point","text"))+labs(title = "PCA (transposée) : projection des échantillons")) #
}




                       #############################
          ############## Fonction PLS-DA Millesime #############
                       #############################



PLSDA_Millesime <- function(){
  #PACKAGES##########
  library(mixOmics) #
  #DONNEES###################################################
  data_filt  <- Pretraitement(FALSE)                        #
  X <- as.matrix(data_filt[ , -c(1:3)])                     #
  Y_ech <- as.factor(sub("\\(.*", "", data_filt$Millesime)) #
  Y_pos <- as.factor(data_filt$Position)                    #
  #ANALYSE#########################################################################################################################
  plsda_model  <- plsda(X, Y_ech, ncomp = 2,scale=TRUE)                                                                           #
  plotIndiv(plsda_model,comp = c(1, 2),group = Y_ech,ind.names = Y_pos,legend = TRUE,title = "PLS-DA : point labels = Position")  #
  plotVar(plsda_model, comp = c(1, 2), title = "Variables sélectionnées (comp1)")                                                 #
  score_vip <- vip(plsda_model)                                                                                                   #
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])                                                    #
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1), ]                                                                    #
  return(vip_df_sorted$mz[1:20])                                                                                                  #
}

                      ##############################
         ############## Fonction sPLS-DA Millesime #############
                      ##############################


SPLSDA_Millesime <- function(bool=TRUE){
  #PACKAGE############
  library(mixOmics)  #
  library(factoextra)#
  #DONNEES###################################################
  data_filt <- Pretraitement(FALSE)                         #
  X <- as.matrix(data_filt[,-c(1:3)])                       #
  Y_ech <- as.factor(sub("\\(.*", "", data_filt$Millesime)) #
  Y_pos <- as.factor(data_filt$Position)                    #
  #ANALYSE####################################################################################################################################################
  splsda_model <- splsda(X, Y_ech, ncomp = 2,keepX = c(20, 20),scale=TRUE)                                                                                   #
  if(bool==TRUE){                                                                                                                                            #
    mixOmics::plotIndiv(splsda_model,comp = c(1, 2),group = Y_ech,ind.names = Y_pos,legend = TRUE,title = "sPLS-DA : point labels = Position",ellipse=FALSE) #
    print(plot_vars(splsda_model))                                                                                                                           #
  }
  score_vip <- vip(splsda_model)                                                                                                  
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])                                                    
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1),]                                                                    
  return(vip_df_sorted$mz[1:20])                                                                                                  
}

                       #####################
          ############## Fonction Plot Var #############
                       #####################


plot_vars <- function(model) {
  #PACKAGES##########
  library(ggplot2)  #
  library(ggrepel)  #
  library(mixOmics) #
  ######################################################################
  # Récupération des coordonnées
  resVar <- plotVar(model, comp = c(1,2), plot = FALSE)
  df <- data.frame(comp1 = resVar$x,comp2 = resVar$y,var   = resVar$names)
  ggplot(df, aes(x = comp1, y = comp2)) +
  geom_path(data = data.frame(theta = seq(0, 2*pi, length.out = 200)),aes(x = cos(theta), y = sin(theta)),inherit.aes = FALSE,linetype = "dashed",colour = "grey80") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 1) +
  geom_text_repel(aes(label = var),size = 2,max.overlaps = Inf,box.padding = 0.5,point.padding = 0.3,segment.size = 0.2,segment.color = "grey50") +
  coord_equal() + labs(title = "Variables sélectionnées (comp1 vs comp2)",x = paste("Composante", c(1)),y = paste("Composante", c(2))) +
  theme_minimal()
}




                            ####################
               ############## Fonction sPLS-DA #############
                            ####################




FonctionSplsDa <- function(bool = TRUE) {
  #PACKAGE ##########
  library(mixOmics) #
  #DONNEES ###########################
  data_filt <- Pretraitement(FALSE)  #
  X <- data_filt[, -c(1:3)]          #
  Y <- as.factor(data_filt$Position) #
  #ANALYSE ########################################################################
  splsda_model <- mixOmics::splsda(X, Y, ncomp = 2, keepX = c(20, 20),scale=TRUE) #
  score_vip <- vip(splsda_model)                                                  #
  vip_df_sorted <- data.frame(mz = rownames(score_vip), VIP1 = score_vip[, 1])    #
  vip_df_sorted <- vip_df_sorted[order(-vip_df_sorted$VIP1),]                     #
  
  if (bool) {
    mixOmics::plotIndiv(splsda_model,comp = c(1,2),group = Y,title = "sPLS-DA (Position)",legend = TRUE,ellipse = TRUE)
    print(plot_vars(splsda_model))
    top_vip_n <- vip_df_sorted[1:20,]
    print(top_vip_n)
  }

  return(vip_df_sorted$mz[1:20])
}


                 ####################
          ######## Fonction ACP Vin #######
                 ####################



PCA_VIN <- function(){
  #PACKAGE ############
  library(readxl)     #
  library(FactoMineR) # 
  library(factoextra) #
  #DONNEES #######################################
  data <- readxl::read_excel("matrice_vin.xlsx") #
  X <- t(data[,-c(1)])                           #
  Y <- as.factor(colnames(data[,-c(1)]))         #
  #ANALYSE ######################################################################################################################################
  res_pca <- FactoMineR::PCA(X,scale.unit = FALSE,graph = FALSE)                                                                                #
  print(fviz_pca_ind(res_pca,habillage = Y, label = "none",mean.point = FALSE,addEllipses = FALSE) +labs(title ="Projection des échantillons")) #
  contrib_total <- rowSums(get_pca_var(res_pca)$contrib[,1:2])                                                                                  #
  top10_total <- names(sort(contrib_total, decreasing = TRUE))[1:10]                                                                            #
  fviz_pca_var(res_pca, select.var = list(name = top10_total),col.var = "contrib")+labs(title ="TOP 10 contribution des variables")             #
}

            ################
     ######## Fonction ACP #######
            ################

PCA <- function(){
  #PACKAGE ############
  library(FactoMineR) # 
  library(factoextra) #
  #DONNEES ######################################################
  data <- Pretraitement(FALSE)                                  #
  data$Millesime <- as.factor(sub("\\(.*", "", data$Millesime)) #
  X <- data[, -c(1:4)]                                          #
  Y <- as.factor(data$Position)                                 #
  Y_ech <- as.factor(data$Millesime)                            #
  #ANALYSE #############################################################################################################################################################
  res_pca <- FactoMineR::PCA(X,scale.unit = TRUE,graph = FALSE)
  print(fviz_pca_ind(res_pca,habillage = Y, label = "none",mean.point = FALSE,addEllipses = FALSE,pointshape  = 19) +labs(title ="Projection des échantillons"))
  print(fviz_pca_ind(res_pca,habillage = Y_ech, label = "none",mean.point = FALSE,addEllipses =FALSE,pointshape  = 19) +labs(title ="Projection des échantillons"))
  
  contrib_total <- rowSums(get_pca_var(res_pca)$contrib[,1:2])
  top10_total <- names(sort(contrib_total, decreasing = TRUE))[1:10]
  fviz_pca_var(res_pca, select.var = list(name = top10_total),col.var = "contrib")+labs(title ="TOP 10 contribution des variables")
}


               #######
    ############ MDS ############
               #######

mds <- function(){
  #PACKAGE #########
  library(ggplot2) #
  #DONNEES ###################
  df <- Pretraitement(FALSE) #
  #ANALYSE #############################################################################################
  mz_cols <- names(df)[!(names(df) %in% c("Millesime","Masse","Position"))]                            #
  cor_mat  <- cor(t(df[, mz_cols]), use = "pairwise.complete.obs")                                     #
  dist_mat <- as.dist(1 - cor_mat)                                                                     #
  mds_res <- cmdscale(dist_mat, k = 2, eig = TRUE)                                                     #
  mds_df <- data.frame(MDS1=mds_res$points[, 1],MDS2=mds_res$points[, 2],Position=factor(df$Position)) #
  #AFFICHAGE ############################################################################################## 
  ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Position)) +                  
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title="MDS des profils m/z",x="Dimension 1",y="Dimension 2",color="Position",shape="Millésime")+
    theme(legend.position = "right")
}

          #####
     ######SNE######
          #####

sne <- function(){
  #PACKAGE #########
  library(ggplot2) #
  library(Rtsne)   #
  #DONNEES #####################
  data_tsne <- Pretraitement() #
  #ANALYSE ###############################################################################
  tsne_res <- Rtsne(data_tsne, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 500) #
  tsne_df <- data.frame(tsne_res$Y)                                                      #
  tsne_df$echantillon <- data_tsne$Millesime                                             #
  tsne_df$Position <- data_tsne$Position                                                 #
  ggplot(tsne_df, aes(x = X1, y = X2, color = echantillon, shape = as.factor(Position))) + geom_point(size = 2) + labs(title = "t-SNE des échantillons", x = "t-SNE 1", y = "t-SNE 2")
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
AffichageBoiteMoustacheMillesime (PLSDA_Millesime())

                   #########################
         ########### sPLS-DA PAR MILLESIME ###########
                   #########################

SPLSDA_Millesime(TRUE)
AffichageBoiteMoustacheMillesime (SPLSDA_Millesime())

                     #######
              ######## PCA #######
                     #######

PCA_VIN()
PCA()


                 #######
      ############ MDS ############
                 #######

mds()

                 #####
          ########SNE########
                 #####

sne()
