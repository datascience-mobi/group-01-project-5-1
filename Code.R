# Daten auslesen und Matrizen erstellen, mit denen gearbeitet wird
`ALLBcell` <- readRDS("ALL-Bcell_list.RDS.gz")
ALLpromotor <- ALLBcell[["promoters"]] 
ALLpromotor <- ALLpromotor[,c(-4,-5,-6,-8,-9,-10)] # Spalten "Strand", "symbol", "entrezID", "GC", "C", "G" gelöscht
ALLpromotorCov <- ALLpromotor[,c(15:24)]  # Nur Coverage Daten
ALLCovMeans <- rowMeans(ALLpromotorCov)   # Liste mit Mittelwert der Coverage pro Genabschnitt
ALLCovMeansLog <- log(ALLCovMeans)        # Bildung Logarithmus für Plot

# Plotten + Speichern der Daten mit Histogramm
png(filename = "E:/Studium HD/4. FS/Bioinfo05/Hist_noLOG.png", width = 1000, height = 500)
hist(ALLCovMeans, xlab = "Means of coverage values", ylab = , xlim = c(0,10000), ylim = c(0,4000), breaks = 1000, lwd = 1) # Nicht logarithmisch
dev.off()

png(filename = "E:/Studium HD/4. FS/Bioinfo05/Hist_LOG.png", width = 1000, height = 500)
hist(ALLCovMeansLog, xlab = "Logarithmic means of coverage values", ylab = , xlim = c(0,13), ylim = c(0,250), breaks = 1000, lwd = 1)  # Logarithmisch
dev.off()

# Threshold bestimmen
quantile(ALLCovMeans, probs = c(0.05, 0.95),na.rm = TRUE) # Gibt für 5% 139, für 95% 5957 -> 5% zu hohe Grenze 
ALLpromotorCov1 <- ALLpromotorCov 
ALLpromotorCov1[ALLpromotorCov > 7753] <- NA # Alle Werte über 98,5% = NA
ALLpromotorCov1[ALLpromotorCov < 28] <- NA  # Alle Werte unter 1,5% = NA

# Isolieren der relevanten Gene (CDH1, CDKN2A, CDKN2B, CDKN1C, KLK10, DKK3, CDH13, PYCARD, DAPK1, PRKN, PTEN, p73, APAF1, LATS1)
ALLCovGen <- ALLpromotorCov[c("ENSG00000039068","ENSG00000078900","ENSG00000147889","ENSG00000147883","ENSG00000129757","ENSG00000129451","ENSG00000050165","ENSG00000140945","ENSG00000103490","ENSG00000120868","ENSG00000196730","ENSG00000185345","ENSG00000131023","ENSG00000171862"),]
ALLCov1Gen <- ALLpromotorCov1[c("ENSG00000039068","ENSG00000078900","ENSG00000147889","ENSG00000147883","ENSG00000129757","ENSG00000129451","ENSG00000050165","ENSG00000140945","ENSG00000103490","ENSG00000120868","ENSG00000196730","ENSG00000185345","ENSG00000131023","ENSG00000171862"),] 

# Überprüfen wie groß der verworfene Anteil
sum(!is.na(ALLpromotorCov)) # Summe aller Variablen in der Matrix die nicht NAs sind
'rejectTotal' <- (1-(sum(!is.na(ALLpromotorCov1))/sum(!is.na(ALLpromotorCov)))) # Prozentsatz der verworfenen Daten insgesammt
'reject' <- (1-(sum(!is.na(ALLCov1Gen))/sum(!is.na(ALLCovGen)))) # Prozentsatz der verworfenen Daten relevanter Gene
rejectTotal # Ausgabe
reject # Ausgabe

# Bei Grenze von 99% bereits 3 Gene mit komplett NA (p73, APAF1, LATS1)
# Grenze letztendlich bei 98,5% (7753) gesetzt: 18% der relevanten Gene verworfen -> 5 Gene (p73, APAF1, LATS1, p57, PTEN) -> bester Kompromiss
# Insgesamt 6,3% aller Daten verworfen (nur über Threshold)

# Übertragen der Coverage NAs auf die Matrix der Beta-values
ALLpromotorBeta <- ALLpromotor[,c(5:14)]  # Erstellen Matrix mit nur Beta values 
for(i in 1:nrow(ALLpromotorCov1)){                      # Abtasten der Zeilen
  for(j in 1:ncol(ALLpromotorCov1)){                    # Abtasten der Spalten
    if(is.na(ALLpromotorCov1[i,j])){    # Prüfen ob NA
      ALLpromotorBeta[i,j] <- NA # Wenn ja, setze Wert bei Beta Values auf NA
    }
  } 
}

repeat{
    # Eliminieren Gene mit >=4 NAs bei den disease Patienten  
    for (i in 1:nrow(ALLpromotorBeta)) {
      for (j in 1:ncol(ALLpromotorBeta)) {       
        if(sum(is.na(ALLpromotorBeta[i,c(1:5)]))>=4){
          ALLpromotorBeta <- ALLpromotorBeta[-i,]
        }
      }
    }
    dim(ALLpromotorBeta) # dim(end) = 55111 (7,9% weg)

    # for (i in 1:nrow(ALLBetaGen)) {  
    #   for (j in 1:ncol(ALLBetaGen)) {       # Nur relevante Gene (Start Anzahl: 14)
    #     if(sum(is.na(ALLBetaGen[i,c(1:5)]))>=4){
    #         ALLBetaGen <- ALLBetaGen[-i,]
    #     }
    #   }
    # }
    # View(ALLBetaGen) # End Anzahl: 9

    # Eliminieren Gene mit >=4 NAs bei den disease Patienten
    for (i in 1:nrow(ALLpromotorBeta)) {
      for (j in 1:ncol(ALLpromotorBeta)) {       
        if(sum(is.na(ALLpromotorBeta[i,c(6:10)]))>=4){
          ALLpromotorBeta <- ALLpromotorBeta[-i,]
        }
      }
    }
    dim(ALLpromotorBeta) 

    # Insgesamt healthy (7,9%) + disease (1%) = 8,9% (5225) aller Gene verworfen 

    # NAs durch Mean der Zeile ersetzen (nach gesund und krank separiert)
    for (i in 1:nrow(ALLpromotorBeta)) {       # Gesunde Daten
      for (j in 1:5) {
         if(is.na(ALLpromotorBeta[i,j])){
           ALLpromotorBeta[i,j] <- rowMeans(ALLpromotorBeta[i,c(1:5)], na.rm = TRUE) 
         }
       }
    }
    for (i in 1:nrow(ALLpromotorBeta)) {       # Kranke Daten
       for (j in 6:10) {
         if(is.na(ALLpromotorBeta[i,j])){
             ALLpromotorBeta[i,j] <- rowMeans(ALLpromotorBeta[i,c(6:10)], na.rm = TRUE) 
         }
       }
    }
    sum(is.na(ALLpromotorBeta))
  if(sum(is.na(ALLpromotorBeta)) == 0){
    break
  }
}

# 0 durch 0,00001 und 1 durch 0,99999 ersetzen
for (i in 1:nrow(ALLpromotorBeta)) {       
  for (j in 1:ncol(ALLpromotorBeta)) {
    if((ALLpromotorBeta[i,j])==0){
      ALLpromotorBeta[i,j] <- 0.00001 
    }
    if((ALLpromotorBeta[i,j])==1){
      ALLpromotorBeta[i,j] <- 0.99999
    }
  }
}

# Normalisierung
ALLMvalue <- ALLpromotorBeta          
for (i in 1:nrow(ALLMvalue)) {
     for (j in 1:ncol(ALLMvalue)) {
         ALLMvalue[i,j] <- log2(ALLMvalue[i,j]/(1-ALLMvalue[i,j])) 
     }
}
View(ALLMvalue)


## PCA
library(ggplot2)
pca <- prcomp(t(ALLMvalue))  # Performs PCA on data matrix and returns results as an object of class prcomp.

# Plotten der Varianz um Daten zu filtern, mit welchen gearbeitet wird
pVar <- (pca$sdev)^2 # Varianz berechnen
plot(pVar, type = "l",main = "Plot of PCA variance", xlab = "PCs") # Elbow bei PC3 -> Nur PC1-3 aussagakräftig
     
# ggPlot
  # Betrachtung der Patienten
  samples <- c(rep("h",5),rep("d",5)) # Erstellen einer weiteren Spalte im Datensatz für Farbigkeit
  tissue <- c(rep("tonsil",3),rep("bone marrow",7))
  provider <- c(rep("J.I.Martin-Subero",3),rep("J.Wiemels",2),rep("A.Bergmann",5))
  date <- c("Aug2013","März2015","Aug2015","Aug2013","März2015","Mai2016","Mai2016",rep("Juli2016",3))
  
  pcax12 <- data.frame(pca$x[,c(1,2)]) # Matrix zu data.frame 
  ggplot(pcax12, aes(x = PC1, y = PC2, colour = samples)) + geom_point() # Erzeugen ggPlot mit farbigen Punkten für PC12, PC13, PC23

  # Alternativ aber hier keine x-/y-Achsenbeschriftung
  # pcx12 <- data.frame(sample=rownames(p$x), X=p$x[,1], Y=p$x[,2]) 
  # ggplot(pcx12, aes(X, Y, group = sample)) + geom_point(aes(color = samples))
  
  # Problem mit Daten: Factor nicht numeric (manuelle Erstellung pcax12-dataframe)
  PC1 <- as.numeric(pca$x[,1])
  PC2 <- as.numeric(pca$x[,2])
  pcax12 <- cbind(PC1,PC2)
  pcax12 <- data.frame(pcax12)
  ggplot(pcax12, aes(x = PC1, y = PC2, colour = samples)) + geom_point()
  
  # Betrachtung der Batch-Effekt: Disease/Healthy, Tissue, Provider, Date
  batcheffect1 <- ggplot(pcax12, aes(x = PC1, y = PC2, colour = samples, shape = provider, size = tissue)) 
  batcheffect1 + geom_point()
  batcheffect2 <- ggplot(pcax12, aes(x = PC1, y = PC2, colour = date, shape = provider, size = tissue))
  batcheffect2 + geom_point()

# Betrachtung der Gene
pcar12 <- data.frame(pca$rotation[,c(1,2)])
ggplot(pcar12, aes(x = PC1, y = PC2)) + geom_point()
     
  # Filtern der Gene mit größter Varianz
  pcaRotVar <- apply(pca$rotation, 1, var) # Varianz der GenPCA
  pcaRotVarsort <- sort(pcaRotVar, decreasing = TRUE) # Sortierung der Varianz absteigen
  plot(pcaRotVarsort, type = "l", main = "Plot der Genvarianz", xlab = "Genes", ylim = c(0,0.0005), xlim = c(0,2000)) # Schauen wo Elbow (hier bei ca. 700)

  # Filtern der Gene nach Loading
  pcaRotAbs <- abs(pca$rotation[,1]) # Betrag der Loadings berechnen
  pcaRotAbssort <- sort(pcaRotAbs, decreasing = TRUE) # Loadings absteigend sortieren
  plot(pcaRotAbssort, type = "l", main = "Plot der PC1 Loading", xlab = "Genes", ylim = c(0,0.06), xlim = c(0,1000)) # Elbow bei unter 50 Genen (seltsam)

  # Genposition nach Varianz ordnen (damit es auf die Methylierungsdaten übertragen werden kann)
  counts <- c(1:54609) # Vektor für Genposition
  pcaRotVar <- cbind(pcaRotVar,counts) # Einbinden des Vektor
  pcaRotVar <- order(pcaRotVar[,1], decreasing = TRUE) # Ordnen der Genpositionen nach Varianz (wie pcaRotVarsort)
  
  # Löschen der Gene hinter Varianz-elbow  
  pcaRotVar <- pcaRotVar[-c(701:54609)] # Alle Gene hinter Elbow (kleinere Varianz als 0.0001388)
  ALLMvalueRemain <- ALLMvalue[pcaRotVar,] # Extrahieren der Mvalues nach den verbliebenen Genpositionen
  ALLMvalueRemain <- rbind(ALLMvalueRemain,ALLMGen) # Miteinbinden der relevanten Gene (sind zuvor nicht im Dataframe gewesen)
      # ALLMvalueRemain: Enthält alle Gene/Patienten nach PCA


# Statistical Tests
pca_remain <- prcomp(t(ALLMvalueRemain)) 
  # Wilcoxon Rank Test
  wilcoxon_p1 <- wilcox.test(pca_remain$x[,1] ~ sample_annotation$TISSUE_TYPE) # Tissuetype von PC1-3 
  wilcoxon_p2 <- wilcox.test(pca_remain$x[,2] ~ sample_annotation$TISSUE_TYPE)
  wilcoxon_p3 <- wilcox.test(pca_remain$x[,3] ~ sample_annotation$TISSUE_TYPE)

  wilcoxon_p4 <- wilcox.test(pca_remain$x[,1] ~ sample_annotation$DISEASE) # Disease PC1-3 
  wilcoxon_p5 <- wilcox.test(pca_remain$x[,2] ~ sample_annotation$DISEASE)
  wilcoxon_p6 <- wilcox.test(pca_remain$x[,3] ~ sample_annotation$DISEASE)

  # Kruskal Wallis Test
  sample_annotation$BIOMATERIAL_PROVIDER <- as.factor(sample_annotation$BIOMATERIAL_PROVIDER) # Sonst Error Gruppenlevel muss endlich sein
  kruskal_p1 <- kruskal.test(pca_remain$x[,1] ~ sample_annotation$BIOMATERIAL_PROVIDER) # Provider PC1-3
  kruskal_p2 <- kruskal.test(pca_remain$x[,2] ~ sample_annotation$BIOMATERIAL_PROVIDER)
  kruskal_p3 <- kruskal.test(pca_remain$x[,3] ~ sample_annotation$BIOMATERIAL_PROVIDER)

  sample_annotation$FIRST_SUBMISSION_DATE <- as.factor(sample_annotation$FIRST_SUBMISSION_DATE)
  kruskal_p4 <- kruskal.test(pca_remain$x[,1] ~ sample_annotation$FIRST_SUBMISSION_DATE) # Date PC1-3 
  kruskal_p5 <- kruskal.test(pca_remain$x[,2] ~ sample_annotation$FIRST_SUBMISSION_DATE)
  kruskal_p6 <- kruskal.test(pca_remain$x[,3] ~ sample_annotation$FIRST_SUBMISSION_DATE)


# Visualisation
  # Heatmap
  Tissue_Type <- c(0.03,0.017,0.83)
  Disease <- c(0.008,0.69,1)
  Provider <- c(0.097,0.097,0.458)
  Submission_Date <- c(0.086,0.617,0.726)
  P_values <- rbind(Tissue_Type, Disease, Provider, Submission_Date) # Erstellen matrix mit p-values
  heatmap(P_values, main = "Heatmap P-values", Colv = NA, Rowv = NA, col = cm.colors(256)) # Heatmap der p-values
  # Levelplot
  library(lattice)
  levelplot(t(P_values), xlab = "PCs", ylab = "Batches", main = "Levelplot P-values") # Levelplot der p-values
  levelplot(t(ALLMvalueRemain[,c(1:3)]), xlab = "PCs", ylab = "Patients", main = "Levelplot PCA") # Levelplot der PCA$x-Werte
  pca_remain_x_abs <- abs(pca_remain$x[,c(1:3)])
  levelplot(t(pca_remain_x_abs), xlab = "PCs", ylab = "Patients", main = "Plot PCA Betrag") # Levelplot des Betrags der PCA$x-Werte
  
# Loading plotten 
pca_rot_PC1 <- pca_remain$rotation[,1] # Verwenden der PC1
pca_rot_PC1_Abs <- abs(pca_rot_PC1) # Loading Betrag berechnen
pca_rot_PC1_sort <- sort(pca_rot_PC1_Abs, decreasing = TRUE) # Loading sortieren
plot(pca_rot_PC1_sort, type = "l",main = "Loading plot PC1", ylab = "Loading", xlab = "Genes")

# K-Means (Silhouette-Plot)
M_km <- kmeans(t(x = ALLMvalueRemain), centers = 2, nstart = 10)
M_km
D <- dist(t(ALLMvalueRemain)) # Distanzmatrix
library(cluster)
silh <- silhouette(M_km$cluster,D) # Silhouettenplot der Patienten
plot(silh, ylab = "Patient") # 2 Cluster erkennbar


# T-test
  Mvalues_noTO <- ALLMvalueRemain[,c(4:10)] # Tonsil-Patienten herausgenommen
  Mvalues_healthy <- Mvalues_noTO[,c(1,2)] # Matrix splitten in healthy und disease
  Mvalues_disease <- Mvalues_noTO[,c(3:7)]

t_test_pvalues <- c(1:709) # Erstellung Liste für p-values

for (i in 1:length(t_test)) {
    t_test[i] <- t.test(Mvalues_healthy[i,],Mvalues_disease[i,])$p.value
}
  # Adjust p-value multiple comparison
  p_correction <- p.adjust(t_test, method = "holm", n = length(t_test))
  p_correction <- order(p_correction) # Korrigierte p-values absteigend ordnen -> Position wird ausgegeben
  t_test <- order(t_test) # Absteigend ordnen -> Position wird angegeben (Gleiche Anordnung wie bei p_correction)

t_test <- t_test[c(1:20)] # Erste 20 Werte behalten
ALLMvalueRemain20 <- ALLMvalueRemain[t_test,] # Dataframe mit verbliebenene M-values
ALLBetaRemain20 <- ALLpromotorBeta[t_test,]   # Dataframe mit verbliebenen Beta-values

# Visualisierung
ALLBetaRemain20 <- as.matrix(ALLBetaRemain20)
levelplot(ALLBetaRemain20, xlab = "Genes", ylab = "Patients", main = "Levelplot Beta-values of remaining Genes")

for (i in 1:20) {
    for(j in 1:10){
        ALLMvalueRemain20[i,j] <- as.numeric(ALLMvalueRemain20[i,j]) # Nur für den Fall, dass Variablen nicht numerisch
    }
}
ALLMvalueRemain20 <- as.matrix(ALLMvalueRemain20)
levelplot(ALLMvalueRemain20, xlab = "Genes", ylab = "Patients", main = "Levelplot M-values of remaining Genes")
