# Daten auslesen und Matrizen erstellen, mit denen gearbeitet wird
`ALLBcell` <- readRDS("ALL-Bcell_list.RDS.gz")
ALLpromotor <- ALLBcell[["promoters"]] 
ALLpromotor <- ALLpromotor[,c(-4,-5,-6,-8,-9,-10)] # Spalten "Strand", "symbol", "entrezID", "GC", "C", "G" gelöscht
ALLpromotorCov <- ALLpromotor[,c(15:24)]  # Nur Coverage Daten
ALLCovMeans <- rowMeans(ALLpromotorCov)   # Liste mit Mittelwert der Coverage pro Genabschnitt
ALLCovMeansLog <- log(ALLCovMeans)        # Bildung Logarithmus für Plot

# Plotten + Speichern der Daten mit Histogramm
png(filename = "E:/Studium HD/4. FS/Bioinfo05/Hist_noLOG.png", width = 1000, height = 500)
hist(ALLCovMeans, xlab = , ylab = , xlim = c(0,10000), ylim = c(0,4000), breaks = 1000, lwd = 1) # Nicht logarithmisch
dev.off()

png(filename = "E:/Studium HD/4. FS/Bioinfo05/Hist_LOG.png", width = 1000, height = 500)
hist(ALLCovMeansLog, xlab = , ylab = , xlim = c(0,13), ylim = c(0,250), breaks = 1000, lwd = 1)  # Logarithmisch
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
           ALLpromotorBeta1[i,j] <- rowMeans(ALLpromotorBeta[i,c(1:5)], na.rm = TRUE) 
         }
       }
    }
    for (i in 1:nrow(ALLpromotorBeta)) {       # Kranke Daten
       for (j in 6:10) {
         if(is.na(ALLpromotorBeta[i,j])){
             ALLpromotorBeta1[i,j] <- rowMeans(ALLpromotorBeta[i,c(6:10)], na.rm = TRUE) 
         }
       }
    }
    sum(is.na(ALLpromotorBeta))
  if(sum(is.na(ALLpromotorBeta)) == 0{
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
p <- prcomp(t(ALLMvalue))  # Performs PCA on data matrix and returns results as an object of class prcomp.

# Plotten der Varianz um Daten zu filtern, mit welchen gearbeitet wird
pVar <- (p$sdev)^2 # Varianz berechnen
plot(pVar, type = "l",main = "Plot of PCA variance", xlab = "PCs") # Elbow bei PC3 -> Nur PC1-3 aussagakräftig
     
# ggPlot
  # Bezogen auf Patienten
samples <- c("h","h","h","h","h","d","d","d","d","d") # Erstellen einer weiteren Spalte im Datensatz für Farbigkeit
p$x <- cbind(p$x, samples) # Einbinden in Matrix
      
pcx12 <- data.frame(p$x[,c(1,2)]) # Matrix zu data.frame 
ggplot(pcx12, aes(x = PC1, y = PC2, colour = samples)) + geom_point() # Erzeugen ggPlot mit farbigen Punkten für PC12, PC13, PC23

# Alternativ aber hier keine x-/y-Achsenbeschriftung
  # pcx12 <- data.frame(sample=rownames(p$x), X=p$x[,1], Y=p$x[,2]) 
  # ggplot(pcx12, aes(X, Y, group = sample)) + geom_point(aes(color = samples))
  
# Datenreduktion anhand der ggPlots
  # Patient 1-3 removed, da anderes Gewebe als Disease
