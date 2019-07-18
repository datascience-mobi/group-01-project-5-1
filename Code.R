# Daten auslesen und Matrizen erstellen, mit denen gearbeitet wird
`ALLBcell` <- readRDS("ALL-Bcell_list.RDS.gz")
ALLpromotor <- ALLBcell[["promoters"]] 
ALLpromotor <- ALLpromotor[,c(-4,-5,-6,-8,-9,-10)] # Spalten "Strand", "symbol", "entrezID", "GC", "C", "G" gelöscht
ALLpromotorCov <- ALLpromotor[,c(15:24)]  # Nur Coverage Daten
ALLCovMeans <- rowMeans(ALLpromotorCov)   # Liste mit Mittelwert der Coverage pro Genabschnitt
ALLCovMeansLog <- log(ALLCovMeans)        # Bildung Logarithmus für Plot

# Plotten der Daten mit Histogramm
hist(ALLCovMeans, xlab = "Means of coverage values", ylab = , xlim = c(0,10000), ylim = c(0,4000), main = "Histogram of the coverage means", breaks = 1000, lwd = 1) # Nicht logarithmisch
hist(ALLCovMeansLog, xlab = "Logarithmic means of coverage values", ylab = , xlim = c(0,13), ylim = c(0,250), main = "Histogram of the logarithmic coverage means",breaks = 1000, lwd = 1)  # Logarithmisch

# Threshold bestimmen
quantile(ALLCovMeans, probs = c(0.05, 0.95),na.rm = TRUE) # Gibt für 5% 139, für 95% 5957 -> 5% zu hohe Grenze 
ALLpromotorCov1 <- ALLpromotorCov 
ALLpromotorCov1[ALLpromotorCov > 7753] <- NA # Alle Werte über 98,5% = NA
ALLpromotorCov1[ALLpromotorCov < 28] <- NA  # Alle Werte unter 1,5% = NA

# Isolieren der selected Genes (CDH1, CDKN2A, CDKN2B, CDKN1C, KLK10, DKK3, CDH13, PYCARD, DAPK1, PRKN, PTEN, p73, APAF1, LATS1)
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
  # Eliminieren Gene mit >=4 NAs bei den healthy Patienten  
  for (i in 1:nrow(ALLpromotorBeta)) {
    for (j in 1:ncol(ALLpromotorBeta)) {       
      if(sum(is.na(ALLpromotorBeta[i,c(1:5)]))>=4){
        ALLpromotorBeta <- ALLpromotorBeta[-i,]
      }
    }
  }
  dim(ALLpromotorBeta) # dim(end) = 55111 (7,9% weg)
  
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
pVarRel <- pVar
pVarRel <- (pVar/sum(pVarRel))*100 # Umrechnen in Prozent
plot(pVarRel, type = "l",main = "Plot of the PCA variance", xlab = "PCs", ylab = "Variance [%]", ylim = c(0,50), xaxt='n')
axis(side = 1, at=1:10) # X-Achsenskalierung ändern

# ggPlot
# Betrachtung der Patienten
samples <- c(rep("h",5),rep("d",5)) # Erstellen einer weiteren Spalte im Datensatz für Farbigkeit
tissue <- c(rep("tonsil",3),rep("bone marrow",7))
provider <- c(rep("J.I.Martin-Subero",3),rep("J.Wiemels",2),rep("A.Bergmann",5))
date <- c("Aug2013","März2015","Aug2015","Aug2013","März2015","Mai2016","Mai2016",rep("Juli2016",3))

pcax12 <- data.frame(pca$x[,c(1,2)]) # Matrix zu data.frame 

# Alternativ aber hier keine x-/y-Achsenbeschriftung
# pcx12 <- data.frame(sample=rownames(p$x), X=p$x[,1], Y=p$x[,2]) 
# ggplot(pcx12, aes(X, Y, group = sample)) + geom_point(aes(color = samples))

# Problem mit Daten: Factor nicht numeric (manuelle Erstellung pcax12-dataframe)
PC1 <- as.numeric(pca$x[,1])
PC2 <- as.numeric(pca$x[,2])
pcax12 <- cbind(PC1,PC2)
pcax12 <- data.frame(pcax12)

# Betrachtung der Batch-Effekt: Disease/Healthy, Tissue, Provider, Date
batcheffect1 <- ggplot(pcax12, aes(x = PC1, y = PC2, colour = samples, shape = provider, size = tissue)) 
batcheffect1 + geom_point()
batcheffect2 <- ggplot(pcax12, aes(x = PC1, y = PC2, colour = date, shape = provider, size = tissue))
batcheffect2 + geom_point()

# Statistical Tests regarding Batch effects
sample_annotation <- read.csv("sample_annotation.csv")
# Wilcoxon Rank Test
wilcoxon_p1 <- wilcox.test(pca$x[,1] ~ sample_annotation$TISSUE_TYPE) # Tissuetype von PC1-3 
wilcoxon_p2 <- wilcox.test(pca$x[,2] ~ sample_annotation$TISSUE_TYPE)
wilcoxon_p3 <- wilcox.test(pca$x[,3] ~ sample_annotation$TISSUE_TYPE)

wilcoxon_p4 <- wilcox.test(pca$x[,1] ~ sample_annotation$DISEASE) # Disease PC1-3 
wilcoxon_p5 <- wilcox.test(pca$x[,2] ~ sample_annotation$DISEASE)
wilcoxon_p6 <- wilcox.test(pca$x[,3] ~ sample_annotation$DISEASE)

# Kruskal Wallis Test
sample_annotation$BIOMATERIAL_PROVIDER <- as.factor(sample_annotation$BIOMATERIAL_PROVIDER) # Sonst Error Gruppenlevel muss endlich sein
kruskal_p1 <- kruskal.test(pca$x[,1] ~ sample_annotation$BIOMATERIAL_PROVIDER) # Provider PC1-3
kruskal_p2 <- kruskal.test(pca$x[,2] ~ sample_annotation$BIOMATERIAL_PROVIDER)
kruskal_p3 <- kruskal.test(pca$x[,3] ~ sample_annotation$BIOMATERIAL_PROVIDER)

sample_annotation$FIRST_SUBMISSION_DATE <- as.factor(sample_annotation$FIRST_SUBMISSION_DATE)
kruskal_p4 <- kruskal.test(pca$x[,1] ~ sample_annotation$FIRST_SUBMISSION_DATE) # Date PC1-3 
kruskal_p5 <- kruskal.test(pca$x[,2] ~ sample_annotation$FIRST_SUBMISSION_DATE)
kruskal_p6 <- kruskal.test(pca$x[,3] ~ sample_annotation$FIRST_SUBMISSION_DATE)


# Visualisation
# Heatmap
Tissue_Type <- c(wilcoxon_p1$p.value,wilcoxon_p2$p.value,wilcoxon_p3$p.value)
Disease <- c(wilcoxon_p4$p.value,wilcoxon_p5$p.value,wilcoxon_p6$p.value)
Provider <- c(kruskal_p1$p.value,kruskal_p2$p.value,kruskal_p3$p.value)
Submission_Date <- c(kruskal_p4$p.value,kruskal_p5$p.value,kruskal_p6$p.value)
P_values <- rbind(Tissue_Type, Disease, Provider, Submission_Date) # Erstellen matrix mit p-values

# Levelplot
library(lattice)
P_values_batch <- P_values
for (i in 1:nrow(P_values)) {
    for (j in 1:ncol(P_values)) {
        if(P_values[i,j] > 0.05){
            P_values_batch[i,j] <- 1
        }
        else{
            P_values_batch[i,j] <- 0
        }
    }
}
levelplot(t(P_values_batch), xlab = "PCs", ylab = "Batches", main = "P-values for different batch effects") # Levelplot der p-values

# Betrachtung der Gene
pcar12 <- data.frame(pca$rotation[,c(1,2)])
ggplot(pcar12, aes(x = PC1, y = PC2)) + geom_point()

# Filtern der Gene mit größter Varianz
pcaRotVar <- apply(pca$rotation, 1, var) # Varianz der GenPCA
pcaRotVarsort <- sort(pcaRotVar, decreasing = TRUE) # Sortierung der Varianz absteigen
plot(pcaRotVarsort, type = "l", main = "Promoter variance for PC1", xlab = "Promoters", ylab ="PC1 variance of promoters", ylim = c(0,0.0005), xlim = c(0,2000)) # Schauen wo Elbow (hier bei ca. 700)
abline(v=700, col = "red")

# Filtern der Gene nach Loading
pcaRotAbs <- abs(pca$rotation[,1]) # Betrag der Loadings berechnen
pcaRotAbssort <- sort(pcaRotAbs, decreasing = TRUE) # Loadings absteigend sortieren
plot(pcaRotAbssort, type = "l", main = "Loading of PC1", xlab = "Promoters", ylab = "PC1 absolute value of promoter loading", ylim = c(0,0.06), xlim = c(0,1000)) # Elbow bei unter 50 Genen (seltsam)
abline(v=30, col = "red")

# Genposition nach Varianz ordnen (damit es auf die Methylierungsdaten übertragen werden kann)
pcaRotVar <- order(pcaRotVar, decreasing = TRUE) # Ordnen der Genpositionen nach Varianz (wie pcaRotVarsort)

# Löschen der Gene hinter Varianz-elbow  
pcaRotVar <- pcaRotVar[-c(701:length(pca$rotation))] # Alle Gene hinter Elbow (kleinere Varianz als 0.0001388)
ALLMvalueRemain <- ALLMvalue[pcaRotVar,] # Extrahieren der Mvalues nach den verbliebenen Genpositionen
ALLMGen <- ALLMvalue[c("ENSG00000039068","ENSG00000147889","ENSG00000147883","ENSG00000129451","ENSG00000050165","ENSG00000140945","ENSG00000103490","ENSG00000196730","ENSG00000185345"),] # relevanten Gene nach Literatur
ALLMvalueRemain <- rbind(ALLMvalueRemain,ALLMGen) # Miteinbinden der relevanten Gene (sind zuvor nicht im Dataframe gewesen)
# ALLMvalueRemain: Enthält alle Gene/Patienten nach PCA

# K-Means (Silhouette-Plot)
M_km <- kmeans(t(x = ALLMvalueRemain), centers = 2, nstart = 10)
M_km$cluster
plot(M_km$cluster, type = "p", xlab = "Patients", ylab = "Clusters", yaxt = "n", pch = 19, main = "K-mean clusters of the samples/patients")
axis(side = 2, at = 1:2)

D <- dist(t(ALLMvalueRemain)) # Distanzmatrix
library(cluster)
silh <- silhouette(M_km$cluster,D) # Silhouettenplot der Patienten
plot(silh, ylab = "Patient", main = "Silhouette plot of patient clustering") # 2 Cluster erkennbar



# T-test
Mvalues_noTO <- ALLMvalueRemain[,c(4:10)] # Tonsil-Patienten herausgenommen
Mvalues_healthy <- Mvalues_noTO[,c(1,2)] # Matrix splitten in healthy und disease
Mvalues_disease <- Mvalues_noTO[,c(3:7)]

t_test <- c(1:709) # Erstellung Liste für p-values

for (i in 1:length(t_test)) {
  t_test[i] <- t.test(Mvalues_healthy[i,],Mvalues_disease[i,])$p.value
}
# Adjust p-value multiple comparison
p_correction <- p.adjust(t_test, method = "holm", n = length(t_test))
p_correction_order <- order(p_correction) # Korrigierte p-values aufsteigend ordnen -> Position wird ausgegeben
p_correction_sort <- sort(p_correction)# Aufsteigend ordnen -> Wert wird angegeben


t_test20 <- p_correction_order[c(1:20)] # Erste 20 Werte behalten der lorrigierten p-values
ALLMvalueRemain20 <- ALLMvalueRemain[t_test20,] # Dataframe mit verbliebenene M-values (20)
ALLBetaRemain20 <- ALLpromotorBeta[t_test20,]   # Dataframe mit verbliebenen Beta-values (20)

threshold <- sum(p_correction_sort <= 0.05) # Threshold auf p = 0.05 gesetzt der korrigierten p-values
t_test_threshold <- p_correction_order[c(1:threshold)]
ALLMvalueRemain_threshold <- ALLMvalueRemain[t_test_threshold,] # Dataframe mit verbliebenene M-values (12)
ALLBetaRemain_threshold <- ALLpromotorBeta[t_test_threshold,]   # Dataframe mit verbliebenen Beta-values (12)

# Visualisierung (Levelplot)
library(lattice)
# Beta-value 20 verbliebene Genes
ALLBetaRemain20_plot <- as.matrix(ALLBetaRemain20) # Konvertierung zur Matrix
rownames(ALLBetaRemain20_plot) <- c(1:nrow(ALLBetaRemain20)) # Anpassen der Zeilennamen
colnames(ALLBetaRemain20_plot) <- c("Patient1_TO_Healthy","Patient2_TO_Healthy","Patient3_TO_Healthy","Patient4_BM_Healthy","Patient5_BM_Healthy","Patient6_BM_Disease","Patient7_BM_Disease","Patient8_BM_Disease","Patient9_BM_Disease","Patient10_BM_Disease")
levelplot(ALLBetaRemain20_plot, xlab = "Promoters", ylab = "Patients", main = "Levelplot Beta-values of 20 remaining promoters") 

# M-value 20 verbliebene Genes
for (i in 1:20) {
  for(j in 1:10){
    ALLMvalueRemain20[i,j] <- as.numeric(ALLMvalueRemain20[i,j]) # Nur für den Fall, dass Variablen nicht numerisch
  }
}
ALLMvalueRemain20_plot <- as.matrix(ALLMvalueRemain20) 
rownames(ALLMvalueRemain20_plot) <- c(1:nrow(ALLMvalueRemain20))
colnames(ALLMvalueRemain20_plot) <- c("Patient1_TO_Healthy","Patient2_TO_Healthy","Patient3_TO_Healthy","Patient4_BM_Healthy","Patient5_BM_Healthy","Patient6_BM_Disease","Patient7_BM_Disease","Patient8_BM_Disease","Patient9_BM_Disease","Patient10_BM_Disease")
levelplot(ALLMvalueRemain20_plot, xlab = "Promoters", ylab = "Patients", main = "Levelplot M-values of 20 remaining promoters") # Levelplot 20 verbliebenen M-values

# Beta-values threshold p=0.05
ALLBetaRemain_threshold_plot <- as.matrix(ALLBetaRemain_threshold)
rownames(ALLBetaRemain_threshold_plot) <- c(1:nrow(ALLBetaRemain_threshold))
colnames(ALLBetaRemain_threshold_plot) <- c("Patient1_TO_Healthy","Patient2_TO_Healthy","Patient3_TO_Healthy","Patient4_BM_Healthy","Patient5_BM_Healthy","Patient6_BM_Disease","Patient7_BM_Disease","Patient8_BM_Disease","Patient9_BM_Disease","Patient10_BM_Disease")
levelplot(ALLBetaRemain_threshold_plot, xlab = "Promoters", ylab = "Patients", main = "Levelplot Beta-values of remaining promoters") # Levelplot verbliebener beta-values hinter 0.05 threshold

# M-values threshold p=0.05
ALLMvalueRemain_threshold_plot <- as.matrix(ALLMvalueRemain_threshold)
rownames(ALLMvalueRemain_threshold_plot) <- c(1:nrow(ALLMvalueRemain_threshold))
colnames(ALLMvalueRemain_threshold_plot) <- c("Patient1_TO_Healthy","Patient2_TO_Healthy","Patient3_TO_Healthy","Patient4_BM_Healthy","Patient5_BM_Healthy","Patient6_BM_Disease","Patient7_BM_Disease","Patient8_BM_Disease","Patient9_BM_Disease","Patient10_BM_Disease")
levelplot(ALLMvalueRemain_threshold_plot, xlab = "Promoters", ylab = "Patients", main = "Levelplot M-values of remaining promoters") # Levelplot verbliebener M-values hinter 0.05 threshold

# Ausgewählte Gene nur zum Vergleich
ALLMGen_plot <- ALLMGen
rownames(ALLMGen_plot) <- c("CDH1","p16","p15","NES-1","DKK-3","CDH13","TMS-1","DAPK-1","PRKN")
colnames(ALLMGen_plot) <- c("Patient1_TO_Healthy","Patient2_TO_Healthy","Patient3_TO_Healthy","Patient4_BM_Healthy","Patient5_BM_Healthy","Patient6_BM_Disease","Patient7_BM_Disease","Patient8_BM_Disease","Patient9_BM_Disease","Patient10_BM_Disease")
levelplot(as.matrix(ALLMGen_plot), xlab = "Genes/Promoters", ylab = "Patients", main = "Levelplot M-values of selected genes/promoters") # Levelplot relevanter Gene

# Logistical regression (threshold bei 99 Genes)
Lg_Mvalues <- data.frame(t(ALLMvalueRemain_threshold)) # Matrix transponieren damit Patienten auf y-Achse
Tumor <- factor(c(rep("0",5),rep("1",5)))
Lg_Mvalues <- cbind(Lg_Mvalues,Tumor) # Einbinden einer Spalte mit Bezeichnung ob Tumor ja/nein

#LGpred <- as.matrix(c(1:10))  # Erstellung Matrix für Ergebnisse
#for (i in 1:nrow(ALLMvalueRemain_threshold)) {
#  glm99 <- glm(Tumor ~ Lg_Mvalues[,i], family = "binomial", data
#               = Lg_Mvalues) # Logistical regression einzelner Gene
#  pred <- as.matrix(predict(glm99, type = "response")) # Anwendung zur Vorhersage
#  LGpred <- cbind(LGpred,pred) # Matrix mit allen Ergebnissen der prediction
#}
#View(LGpred)
#LGpred <- LGpred[,-1]
#
#rownames(LGpred) <- c("Patient1_TO_Healthy","Patient2_TO_Healthy","Patient3_TO_Healthy","Patient4_BM_Healthy","Patient5_BM_Healthy","Patient6_BM_Disease","Patient7_BM_Disease","Patient8_BM_Disease","Patient9_BM_Disease","Patient10_BM_Disease")
#levelplot(t(LGpred), xlab = "Promoters", ylab = "Patients", main = "Levelplot of logistical regression prediction") # Visualisierung

Lg_Mvalues2 <- Lg_Mvalues
healthy = which(Lg_Mvalues2$Tumor=="0")
healthy.train = sample(healthy,floor(0.75*length(healthy)))
cancer = which(Lg_Mvalues2$Tumor=="1")
cancer.train = sample(cancer, floor(0.75*length(cancer)))
train = c(healthy.train,cancer.train)
train.set = Lg_Mvalues2[train, 1:10]
test.set = Lg_Mvalues2[-train, 1:10]
Tumor2 <- c(0,0,0,1,1,1)
regression_model2 <- glm(formula = Tumor2 ~ ., family = binomial(link = "logit"), data = train.set)
predict2 <- predict(regression_model2, newdata = test.set, type = "response")

levelplot(t(predict2), xlab = "Promoters", ylab = "Patients", main = "Levelplot of logistical regression prediction")
