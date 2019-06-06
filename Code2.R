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
quantile(ALLMean, probs = c(0.05, 0.95),na.rm = TRUE) # Gibt für 5% 139, für 95% 5957 -> 5% zu hohe Grenze 
ALLpromotorCov1[ALLpromotorCov >= 5957] <- NA # Alle Werte über 95% = NA
ALLpromotorCov1[ALLpromotorCov1 <= 28] <- NA  # Alle Werte unter 5% = NA

# Überprüfen wie groß der verworfene Anteil
sum(!is.na(ALLpromotorCov)) # Summe aller Variablen in der Matrix die nicht NAs sind
'verworfen' <- (1-(sum(!is.na(ALLpromotorCov1))/sum(!is.na(ALLpromotorCov)))) # Prozentsatz der Verworfenen Daten

# Übertragen der Coverage NAs auf die Matrix der Beta-values
ALLpromotor[,c(5:14)] -> ALLpromotorBeta # Erstellen Matrix mit nru Beta values 
for(i in 1:59812){                      # Abtasten der Zeilen
     for(j in 1:10){                    # Abtasten der Spalten
         if(is.na(ALLpromotorCov1[i,j])){    # Prüfen ob NA
             ALLpromotorBeta[i,j] <- NA # Wenn ja, setze Wert bei Beta Values auf NA
         }
     } 
}