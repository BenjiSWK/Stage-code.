library(readxl)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)

####################################
### OUVERTURE DES JEU DE DONNEES ###
####################################
#
Data_1 <- read_excel("Extraction méthanol verticale Comtes Lafon.xlsx")
Data_2 <- read_excel("All Matrice sorted age.xls")

########################################
### NETTOYAGE PREMIER JEU DE DONNEES ###
########################################
#

Data_1_simplifiee = Data_1[c(-1,-2,-3,-4,-5,-6),c(-3,-4)]
names(Data_1_simplifiee)[2] <- "Position"
names(Data_1_simplifiee)[3] <- "Masse"
names(Data_1_simplifiee)[1] <- "Millesime"
Data_1_simplifiee <- fill(Data_1_simplifiee,Millesime, .direction = "down")
Data_1_simplifiee <- Data_1_simplifiee %>%
  mutate(
    Millesime = str_replace_all(Millesime, " ", ""),
    echantillon_position = paste0(
      # Début de la clé
      ifelse(
        str_detect(Position, "M"),
        paste0(Millesime, "-", Position),
        paste0(Millesime, "-", Position, "M")
      ),
      # Ajout du suffixe
      "_FORMULAE.dat"
    )
  )     

#On supprime les échantillons E1 et E2 du millesime 1969 (voir le rapport)
Data_1_simplifiee$Position <- gsub("MI", "", Data_1_simplifiee$Position)
lignes_a_exclure <- grepl("ME1|ME2", Data_1_simplifiee$Position)
Data_1_simplifiee <- Data_1_simplifiee[!lignes_a_exclure, ]
#On vient supprimer les bouteilles bandol et cork
lignes_a_exclure <- grepl("cork|bandol", Data_1_simplifiee$Millesime, ignore.case = TRUE)
Data_1_simplifiee <- Data_1_simplifiee[!lignes_a_exclure, ]

#########################################
### NETTOYAGE DEUXIEME JEU DE DONNEES ###
#########################################
#
Data_2_simplifiee <- as.data.frame(t(Data_2[,-1]))
Data_2_simplifiee <- cbind(rownames(Data_2_simplifiee), Data_2_simplifiee)
names(Data_2_simplifiee)[1] <- "echantillon_position"
rownames(Data_2_simplifiee) <- NULL

################
### JOINTURE ###
################
#
data_combine <- Data_2_simplifiee %>% right_join(Data_1_simplifiee, by = "echantillon_position") 
data_combine <- data_combine %>%select(Millesime,Masse,Position, everything())  
data_combine <- data_combine %>% mutate(across(-c(Millesime,Position,echantillon_position), as.numeric))
colnames(data_combine)[-c(1:4)] <- Data_2$`Mass (avg.)`

####################
### EXPORTATION  ###
####################
#
data_combine$echantillon_position <- NULL
library(writexl)
data_n <- data_combine
write_xlsx(data_n, "data_n.xlsx")
