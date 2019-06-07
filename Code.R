# Daten auslesen und Matrizen erstellen, mit denen gearbeitet wird
`ALLBcell` <- readRDS("E:/Studium/4. FS/Bioinfo05/project-05-group-01/ALL-Bcell_list.RDS")
ALLpromotor[,c(15:24)] -> ALLpromotorCov # Nur Coverage Daten
rowMeans(ALLpromotorCov) -> ALLCovMeans  # Liste mit Mittelwert der Coverage pro Genabschnitt
log(ALLCovMeans) -> ALLCovMeansLog       # Bildung Logarithmus für Plot

# Plotten + Speichern der Daten mit Histogramm
png(filename = "E:/Studium HD/4. FS/Bioinfo05/Hist_noLOG.png", width = 1000, height = 500)
hist(ALLCovMeans, xlab = , ylab = , xlim = c(0,10000), ylim = c(0,4000), breaks = 1000, lwd = 1) # Nicht logarithmisch
dev.off()

png(filename = "E:/Studium HD/4. FS/Bioinfo05/Hist_LOG.png", width = 1000, height = 500)
hist(ALLCovMeansLog, xlab = , ylab = , xlim = c(0,13), ylim = c(0,250), breaks = 1000, lwd = 1)  # Logarithmisch
dev.off()

# Threshold bestimmen
quantile(ALLCovMeans, probs = c(0.05, 0.95),na.rm = TRUE) # Gibt für 5% 139, für 95% 5957 -> 5% zu hohe Grenze 
ALLpromotorCov -> ALLpromotorCov1
ALLpromotorCov1[ALLpromotorCov >= 8342] <- NA # Alle Werte über 99% = NA
ALLpromotorCov1[ALLpromotorCov1 <= 28] <- NA  # Alle Werte unter 1,5% = NA

# Übertragen der Coverage NAs auf die Matrix der Beta-values
ALLpromotor[,c(5:14)] -> ALLpromotorBeta # Erstellen Matrix mit nru Beta values 
for(i in 1:59812){                      # Abtasten der Zeilen
     for(j in 1:10){                    # Abtasten der Spalten
         if(is.na(ALLpromotorCov1[i,j])){    # Prüfen ob NA
             ALLpromotorBeta[i,j] <- NA # Wenn ja, setze Wert bei Beta Values auf NA
         }
     } 
}

# Selbst bei Grenze von 99% 3 Gene mit komplett NA (p73, APAF1, LATS1) deshalb diese entfernt

# Isolieren der relevanten Gene (CDH1, CDKN2A, CDKN2B, CDKN1C, KLK10, DKK3, CDH13, PYCARD, DAPK1, PRKN, PTEN)
ALLpromotorCov[c("ENSG00000039068","ENSG00000078900","ENSG00000147889","ENSG00000147883","ENSG00000129757","ENSG00000129451","ENSG00000050165","ENSG00000140945","ENSG00000103490","ENSG00000120868","ENSG00000196730","ENSG00000185345","ENSG00000131023","ENSG00000171862"),] -> ALLCovGen
ALLpromotorCov1[c("ENSG00000039068","ENSG00000078900","ENSG00000147889","ENSG00000147883","ENSG00000129757","ENSG00000129451","ENSG00000050165","ENSG00000140945","ENSG00000103490","ENSG00000120868","ENSG00000196730","ENSG00000185345","ENSG00000131023","ENSG00000171862"),] -> ALLCov1Gen
ALLpromotorBeta[c("ENSG00000039068","ENSG00000078900","ENSG00000147889","ENSG00000147883","ENSG00000129757","ENSG00000129451","ENSG00000050165","ENSG00000140945","ENSG00000103490","ENSG00000120868","ENSG00000196730","ENSG00000185345","ENSG00000131023","ENSG00000171862"),] -> ALLBetaGen

# Überprüfen wie groß der verworfene Anteil
sum(!is.na(ALLpromotorCov)) # Summe aller Variablen in der Matrix die nicht NAs sind
'verworfenGes' <- (1-(sum(!is.na(ALLpromotorCov1))/sum(!is.na(ALLpromotorCov)))) # Prozentsatz der verworfenen Daten insgesammt
'verworfen' <- (1-(sum(!is.na(ALLCov1Gen))/sum(!is.na(ALLCovGen)))) # Prozentsatz der verworfenen Daten relevanter Gene
verworfenGes # Ausgabe
verworfen # Ausgabe
View(ALLBetaGen) # Direkt schauen wie viele Gene/Werte NA

# Grenze letztendlich bei 98,5% (7753) gesetzt, da 10% der Daten verworfen + 2 weitere Gene verworfen (p57, PTEN) -> bester Kompromiss

# Eliminieren Gene mit zu vielen NAs
ALLBetaGen[c(-4,-11),] -> ALLBetaGen # Alternative finden mit Zeilen, welche NA >= 4 rausschmeißen

# Normalisierung
ALLBetaGen -> ALLMGen
for (j in 1:10) {
     for (i in 1:9) {
         log2(ALLBetaGen[i,j]/(1-ALLBetaGen[i,j])) -> ALLMGen[i,j]
     }
}
View(ALLMGen)

# NAs durch Mean der Zeile ersetzen
for (j in 1:10) {
     for (i in 1:9) {
         if(is.na(ALLMGen[i,j])){
             rowMeans(ALLMGen[i,], na.rm = TRUE) -> ALLMGen[i,j]
         }
     }
 }
> View(ALLMGen)