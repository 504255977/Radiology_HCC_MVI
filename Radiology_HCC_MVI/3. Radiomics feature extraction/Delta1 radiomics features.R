library(tidyverse)

##### Delta1 radiomics features ####
data <- read.csv('radiomics features/zhejiang tumor.csv')

shape <- data[,c(1:15)] #shape radiomics features
n <- data[,c(16:107)] #NC original radiomics features 
a <- data[,c(108:199)] #Arterial phase original radiomics features 
v <- data[,c(200:291)] #Portal phase original radiomics features
d <- data[,c(292:383)] #Delayed phase original radiomics features

an <- as.data.frame(data.matrix(a)-data.matrix(n))
vn <- as.data.frame(data.matrix(v)-data.matrix(n))
dn <- as.data.frame(data.matrix(d)-data.matrix(n))
va <- as.data.frame(data.matrix(v)-data.matrix(a))
da <- as.data.frame(data.matrix(d)-data.matrix(a))
dv <- as.data.frame(data.matrix(d)-data.matrix(v))

data_output <- shape%>%
  cbind(an)%>%
  cbind(vn)%>%
  cbind(dn)%>%
  cbind(va)%>%
  cbind(da)%>%
  cbind(dv)

#names(data_output) <- names(XXXXXXXX) #Define names of delta1 radiomics features
write.csv(data_output,file = "radiomics features/zhejiang tumor delta1.csv",row.names = FALSE)
