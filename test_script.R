
filename <- "C:/Users/erik_/Documents/erik documents/Programas/Inferencia_Bayesiana_Proyecto_Final/INFERENCIA.txt"

# Read file
df <- read.table(file = filename, header = TRUE ,sep = ' ',skip = 1 )
head(df)

df["deltacp"]
