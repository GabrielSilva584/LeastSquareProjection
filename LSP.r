# Bibliotecas usadas
# Necessário instalar kmed
# install.packages("kmed")
# Necessário instalar FNN
# install.packages("FNN")
library("MASS")
library("FNN")
library("kmed")

lsp <- function(data){
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
    
    X.medoid <- matrix(data=0, nrow=num_cl, ncol=2)
    X.medoid[1:num_cl,] <- X[cl$medoid,]
    
    return(list(points=X, medoids=X.medoid))
}

# Pegar argumentos de linha de comando
args = commandArgs(trailingOnly=TRUE)

plot.path <- "plot.pdf"

if (length(args)>=1) {
    # Ler o arquivo de entrada data.csv
    data <- read.csv(arg[1], header = TRUE, sep = ",", dec = ".")
    if(length(args)>=2){
        plot.path <- args[2]
    }
}else{
    data("mtcars", "iris")
    data <- mtcars
    
    iris$Species <- NULL
    X2 <- lsp(iris)
}

X <- lsp(data)

# Plot simples
pdf(plot.path)
plot(X$points, col="blue", pch=16)
points(X$medoids,col="red", pch=16)

if(exists("X2")){
    plot(X2$points, col="blue", pch=16)
    points(X2$medoids,col="red", pch=16)
}

dev.off()
