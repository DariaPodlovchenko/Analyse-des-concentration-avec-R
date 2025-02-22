# Importation des données
file_path <- "PM10_NORD_2019-2020.csv"

# Lecture des données depuis un fichier CSV avec un séparateur ';'
data <- read.csv(file_path, header = TRUE, sep = ";")

# Renommage des colonnes
colnames(data) <- c("Department", "StationCode", "Typology", "Pollutant", 
                    "PM10_2019", "PM10_2020", "Longitude", "Latitude")

# Taille des données
dim(data)

# Noms des colonnes
colnames(data)

# Types de données 
sapply(data, class)

# Premiers enregistrements du tableau
head(data, n = 5)

# Vérification des valeurs manquantes
colSums(is.na(data))

# Résumé 
summary(data[, c("PM10_2019", "PM10_2020")])

# Création de boxplots
boxplot(data$PM10_2019, data$PM10_2020,
        names = c("2019", "2020"),  # Légendes des boxplots
        main = "Comparaison des PM10 (2019 vs 2020)",  # Titre
        ylab = "Concentration de PM10 (µg/m3)",  # Légende de l'axe Y
        col = c("skyblue", "salmon"),  # Couleurs pour chaque année
        border = "black")  # Couleur des bordures

# Histogrammes
par(mfrow = c(1, 2))  # Fenêtre graphique divisée en 2 parties (1 ligne, 2 colonnes)

# Histogramme pour 2019
hist(data$PM10_2019,
     main = "Distribution de PM10 (2019)",  # Titre
     xlab = "Concentration de PM10 (µg/m3)",  # Légende de l'axe X
     ylab = "Fréquence",  # Légende de l'axe Y
     col = "skyblue",  # Couleur de l'histogramme
     xlim = c(10, 30),  # Échelle commune pour comparaison
     breaks = seq(10, 30, by = 2))  # Intervalles

# Histogramme pour 2020
hist(data$PM10_2020,
     main = "Distribution de PM10 (2020)",  # Titre
     xlab = "Concentration de PM10 (µg/m3)",  # Légende de l'axe X
     ylab = "Fréquence",  # Légende de l'axe Y
     col = "salmon",  # Couleur de l'histogramme
     xlim = c(10, 30),  # Échelle commune pour comparaison
     breaks = seq(10, 30, by = 2))  # Intervalles


# 3. Test de significativité des différences
# Test de la différence des moyennes (t-test)
t_test <- t.test(data$PM10_2019, data$PM10_2020, paired = TRUE)
cat("Résultat du t-test:\n")
print(t_test)

# Test de la différence des variances (F-test)
var_test <- var.test(data$PM10_2019, data$PM10_2020)
print(var_test)

# 4. Test de corrélation linéaire
cor_test <- cor.test(data$PM10_2019, data$PM10_2020, method = "pearson")
print(cor_test)

# 5. Analyse de l'autocorrélation spatiale
library(vegan)
library(ape)  
coords <- as.matrix(data[, c("Longitude", "Latitude")]) 

# Matrice des distances
DistMat <- as.matrix(dist(coords))  # Distances euclidiennes
DistMat[DistMat == 0] <- Inf  # Évite la division par 0

# Matrice des poids
WeightMat <- 1 / DistMat  # Poids inversement proportionnels aux distances
diag(WeightMat) <- 0  # Pas d'auto-influence

# Calcul de l'indice de Moran
moran_2019 <- Moran.I(data$PM10_2019, WeightMat)
moran_2020 <- Moran.I(data$PM10_2020, WeightMat)
print(moran_2019)
print(moran_2020)

# Test de Mantel
Dist_PM10_2019 <- as.matrix(dist(data$PM10_2019))
Dist_PM10_2020 <- as.matrix(dist(data$PM10_2020))
mantel_2019 <- mantel(DistMat, Dist_PM10_2019)
mantel_2020 <- mantel(DistMat, Dist_PM10_2020)
print(mantel_2019)
print(mantel_2020)

# Calcul de l'indice de Geary
library(pgirmess)  

# Corrélogramme et indice de Geary pour 2019
geary_2019 <- correlog(coords, data$PM10_2019, method = "Geary")
print(geary_2019)

# Corrélogramme et indice de Geary pour 2020
geary_2020 <- correlog(coords, data$PM10_2020, method = "Geary")
print(geary_2020)

# Visualisation des indices de Geary
plot(geary_2019, main = "Indice de Geary pour 2019")
plot(geary_2020, main = "Indice de Geary pour 2020")


# 6. Modélisation de la structure spatiale avec Krigeage et analyse des performances

library(fields)

# Création d'une grille régulière pour les prédictions spatiales
grid_lon <- seq(min(data$Longitude), max(data$Longitude), length.out = 100)
grid_lat <- seq(min(data$Latitude), max(data$Latitude), length.out = 100)
grid <- expand.grid(Lon = grid_lon, Lat = grid_lat)

# Modélisation par krigeage pour 2019
krig_2019 <- Krig(coords, data$PM10_2019)  # Ajustement du modèle
interp_2019 <- predictSurface(krig_2019, list(x = grid_lon, y = grid_lat))  # Prédictions

# Modélisation par krigeage pour 2020
krig_2020 <- Krig(coords, data$PM10_2020)  # Ajustement du modèle
interp_2020 <- predictSurface(krig_2020, list(x = grid_lon, y = grid_lat))  # Prédictions

# Évaluation des modèles avec les données réelles (correlation test)
pred_2019 <- predict(krig_2019, coords)  # Prédictions pour les stations en 2019
pred_2020 <- predict(krig_2020, coords)  # Prédictions pour les stations en 2020

# Calcul des corrélations entre prédictions et données observées
cor_test_2019 <- cor.test(data$PM10_2019, pred_2019)
cor_test_2020 <- cor.test(data$PM10_2020, pred_2020)

# Vérification de la normalité des résidus avec le test de Shapiro-Wilk
residuals_2019 <- data$PM10_2019 - pred_2019
residuals_2020 <- data$PM10_2020 - pred_2020

shapiro_test_2019 <- shapiro.test(residuals_2019)
shapiro_test_2020 <- shapiro.test(residuals_2020)

# Affichage des résultats
print(cor_test_2019)
print(cor_test_2020)
print(shapiro_test_2019)
print(shapiro_test_2020)

# Visualisation des résidus
par(mfrow = c(1, 2))  # Diviser l'affichage en 2 colonnes
hist(residuals_2019, main = "Résidus 2019", xlab = "Résidus", col = "skyblue", breaks = 10)
hist(residuals_2020, main = "Résidus 2020", xlab = "Résidus", col = "salmon", breaks = 10)

# Visualisation des cartes d'interpolation
par(mfrow = c(1, 2))

# Carte pour 2019
image.plot(interp_2019$x, interp_2019$y, interp_2019$z,
           main = "Interpolation PM10 (2019)",
           xlab = "Longitude", ylab = "Latitude", col = topo.colors(100))
points(coords[, 1], coords[, 2], pch = 19, col = "red")
contour(interp_2019$x, interp_2019$y, interp_2019$z, add = TRUE, col = "blue")

# Carte pour 2020
image.plot(interp_2020$x, interp_2020$y, interp_2020$z,
           main = "Interpolation PM10 (2020)",
           xlab = "Longitude", ylab = "Latitude", col = topo.colors(100))
points(coords[, 1], coords[, 2], pch = 19, col = "red")
contour(interp_2020$x, interp_2020$y, interp_2020$z, add = TRUE, col = "blue")



# Carte finale
ModGrid_2019 <- predictSurface(krig_2019, list(x = grid_lon, y = grid_lat))
ModGrid_2019$z[ModGrid_2019$z < 14] <- 14  # Limite basse pour correspondre à zlim
ModGrid_2019$z[ModGrid_2019$z > 22] <- 22  # Limite haute pour correspondre à zlim

# Utilisation des prédictions sur une surface pour 2020
ModGrid_2020 <- predictSurface(krig_2020, list(x = grid_lon, y = grid_lat))
ModGrid_2020$z[ModGrid_2020$z < 14] <- 14  # Limite basse pour correspondre à zlim
ModGrid_2020$z[ModGrid_2020$z > 22] <- 22  # Limite haute pour correspondre à zlim

# Visualisation des deux années avec les mêmes paramètres
par(mfrow = c(1, 2))  # Deux graphiques côte à côte

# Carte pour 2019
surface(ModGrid_2019, xlab = "Longitude", ylab = "Latitude",
        main = "Carte des PM10 (2019)", zlim = c(14, 22), col = topo.colors(100))
map("france", add = TRUE)  # Ajout des limites des départements

# Carte pour 2020
surface(ModGrid_2020, xlab = "Longitude", ylab = "Latitude",
        main = "Carte des PM10 (2020)", zlim = c(14, 22), col = topo.colors(100))
map("france", add = TRUE)  # Ajout des limites des départements
