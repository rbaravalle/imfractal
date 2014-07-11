library("kohonen")

breadtrain <- read.csv("breadtrainS.csv", sep=",", h=F)
breadtest <- read.csv("breadtestS.csv", sep=",", h=F)
nonbreadtrain <- read.csv("nonbreadtrainS.csv", sep=",", h=F)
nonbreadtest <- read.csv("nonbreadtestS.csv", sep=",", h=F)

breads <- rbind(breadtrain,breadtest,nonbreadtrain,nonbreadtest)

labels <- floor(breads[1])
bbread<- breads[-1]

set.seed(1)
     
training <- sample(nrow(bbread), 32)
Xtraining <- scale(bbread[training,])

grid <- somgrid(4, 4, "hexagonal")
som.data <- som(Xtraining, grid)
summary(som.data)

set.seed(4)
png("som_sandbox.png")
kohonen::plot.kohonen(som.data, type = "mapping", label = labels[training,], col = labels[training,], main = "")
dev.off()

q()


breads <- read.csv("som2.txt", sep=" ", h=F)

set.seed(2)
     
training <- sample(nrow(breads), 100)
Xtraining <- scale(breads[training, -4])

grid = somgrid(10, 10, "hexagonal")
som.breads <- som(Xtraining, grid)


png("som.rgb.png")
kohonen::plot.kohonen(som.breads, type = "mapping", label = breads[training,4], col = breads[training,4]+1, main="")
dev.off()
