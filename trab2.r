# Bibliotecas usadas
# Necessário instalar kmed
# install.packages("kmed")
# Necessário instalar FNN
# install.packages("FNN")
# Necessário instalar highcharter (requer gfortran)
# install.packages("highcharter")
library("MASS")
library("FNN")
library("kmed")
library("highcharter")

# Ler o arquivo de entrada data.csv
file = "data.csv"
data <- read.csv(file, header = TRUE, sep = ",", dec = ".")
data.matrix <- as.matrix(data)

# Propriedades dos dados
size <- nrow(data)
num_atrib <- ncol(data)
num_cl = sqrt(size)

# Calcular os clusters por k-medoids
cl <- fastkmed(dist(data), ncluster=num_cl)
data.medoid <- matrix(data=0, nrow=num_cl, ncol=num_atrib)
data.medoid[1:num_cl,] <- data.matrix[cl$medoid,]
data.medoid.sammon <- sammon(dist(data.medoid))

# Calcular as vizinhanças por knn
knn <- get.knn(data, k=num_cl)

# Criar a matriz L
L <- matrix(data=0, nrow=size, ncol=size)

# Preencher a matriz L
for(i in 1:size){
    L[i, i] <- 1
    L[i, knn$nn.index[i,]] <- -1/num_cl
}

# Preencher a matriz C
C <- matrix(data=0, nrow=num_cl, ncol=size)
for(i in 1:num_cl){
    C[i, cl$medoid[i]] <- 1
}

# Criar matriz A
A <- rbind(L,C)
At <- t(A)

# Preencher matriz B
B <- matrix(data=0, nrow=size+num_cl, ncol=2)
B[(size+1):(size+num_cl),] <- data.medoid.sammon$points[1:num_cl,]

# Solucionar o sistema
X <- solve(At %*% A, At %*% B)

# Plot simples
plot(data.medoid.sammon$points, col="red", pch=15)
points(X, col="blue", pch=16)

# Plot interativo
hchart(X, "scatter")
