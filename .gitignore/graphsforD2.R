library(igraph)
library(Matrix)
library(lattice)
library(vcd)

netmat1 <- rbind(c(1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0),
                 c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0),
                 c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0),
                 c(0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
                 c(0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0),
                 c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0),
                 c(0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0),
                 c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0),
                 c(0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0),
                 c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0),
                 c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
                 c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
                 c(1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
rownames(netmat1) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
colnames(netmat1) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
netmat1
class(netmat1)

isSymmetric(netmat1)  #YES

net1<- graph_from_adjacency_matrix( netmat1, mode = "undirected",weighted = NULL, 
  diag = FALSE, add.colnames = NULL, add.rownames = NA) 
net1
class(net1)
laynet<-layout_in_circle(net1)
plot(net1, layout=laynet, vertex.frame.color="royalblue", vertex.color="royalblue",
     vertex.size=2.5, edge.width=0.5)
title(main=c("Νetwork from Pearson Corelation for the 7th window"))
dev.off()


netmat2 <- rbind(c(1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1),
                   c(0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0),
                   c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0),
                   c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0),
                   c(1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1),
                   c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0),
                   c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
                   c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0),
                   c(1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1),
                   c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0),
                   c(0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1),
                   c(1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1))
rownames(netmat2) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
colnames(netmat2) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
netmat2
class(netmat2)

isSymmetric(netmat2)  #YES

net2<- graph_from_adjacency_matrix( netmat2, mode = "undirected",weighted = NULL, 
                                    diag = FALSE, add.colnames = NULL, add.rownames = NA) 
net2
class(net2)
laynet<-layout_in_circle(net2)
plot(net2, layout=laynet, vertex.frame.color="royalblue", vertex.color="royalblue",
     vertex.size=2.5, edge.width=0.5)
title

######################################################
###Normality matrices
netmat3 <- rbind(c( 1.0        ,  0.46666667,  0.86666667,  0.4       ,  0.13333333,
                    0.26666667,  0.26666667,  0.26666667,  0.66666667,  0.33333333,
                    0.26666667,  0.26666667,  0.46666667,  0.26666667,  0.13333333,
                    0.        ,  0.26666667,  0.4       ,  0.13333333,  0.6       ),
                  c( 0.46666667,  1.0        ,  0.86666667,  0.26666667,  0.4       ,
                    0.4       ,  0.13333333,  0.53333333,  0.66666667,  0.6       ,
                    0.13333333,  0.46666667,  0.6       ,  0.13333333,  0.13333333,
                    0.0        ,  0.13333333,  0.4       ,  0.4       ,  0.2       ),
                  c( 0.86666667,  0.86666667,  1.        ,  0.66666667,  0.8       ,
                    0.53333333,  0.4       ,  0.66666667,  0.53333333,  0.46666667,
                    0.4       ,  0.33333333,  0.33333333,  0.26666667,  0.13333333,
                    0.        ,  0.13333333,  0.53333333,  0.66666667,  0.2       ),
                  c(0.4       ,  0.26666667,  0.66666667,  1.        ,  0.8       ,
                    1.        ,  0.53333333,  0.53333333,  0.53333333,  0.53333333,
                    0.53333333,  0.4       ,  0.13333333,  0.4       ,  0.26666667,
                    0.13333333,  0.        ,  0.4       ,  0.66666667,  0.4       ),
                  c(0.13333333,  0.4       ,  0.8       ,  0.8       ,  1.        ,
                    0.53333333,  0.8       ,  0.53333333,  0.13333333,  0.4       ,
                    0.4       ,  0.26666667,  0.13333333,  0.2       ,  0.2       ,
                    0.        ,  0.        ,  0.26666667,  0.53333333,  0.13333333),
                  c(0.26666667,  0.4       ,  0.53333333,  1.        ,  0.53333333,
                    1.        ,  0.4       ,  0.4       ,  0.4       ,  0.4       ,
                    0.53333333,  0.26666667,  0.13333333,  0.26666667,  0.4       ,
                    0.        ,  0.2       ,  0.4       ,  0.6       ,  0.4       ),
                  c(0.26666667,  0.13333333,  0.4       ,  0.53333333,  0.8       ,
                    0.4       ,  1.        ,  0.33333333,  0.26666667,  0.13333333,
                    0.26666667,  0.4       ,  0.13333333,  0.33333333,  0.2       ,
                    0.        ,  0.        ,  0.33333333,  0.46666667,  0.26666667),
                  c(0.26666667,  0.53333333,  0.66666667,  0.53333333,  0.53333333,
                    0.4       ,  0.33333333,  1.        ,  0.4       ,  0.53333333,
                    0.53333333,  0.53333333,  0.13333333,  0.46666667,  0.46666667,
                    0.        ,  0.        ,  0.46666667,  0.26666667,  0.        ),
                  c(0.66666667,  0.66666667,  0.53333333,  0.53333333,  0.13333333,
                    0.4       ,  0.26666667,  0.4       ,  1.        ,  0.53333333,
                    0.53333333,  0.53333333,  0.53333333,  0.4       ,  0.4       ,
                    0.26666667,  0.06666667,  0.46666667,  0.6       ,  0.13333333),
                  c(0.33333333,  0.6       ,  0.46666667,  0.53333333,  0.4       ,
                    0.4       ,  0.13333333,  0.53333333,  0.53333333,  1.        ,
                    0.66666667,  0.86666667,  0.33333333,  0.4       ,  0.26666667,
                    0.        ,  0.13333333,  0.26666667,  0.53333333,  0.33333333),
                  c(0.26666667,  0.13333333,  0.4       ,  0.53333333,  0.4       ,
                    0.53333333,  0.26666667,  0.53333333,  0.53333333,  0.66666667,
                    1.        ,  0.26666667,  0.26666667,  0.4       ,  0.26666667,
                    0.13333333,  0.        ,  0.4       ,  0.4       ,  0.33333333),
                  c(0.26666667,  0.46666667,  0.33333333,  0.4       ,  0.26666667,
                    0.26666667,  0.4       ,  0.53333333,  0.53333333,  0.86666667,
                    0.26666667,  1.        ,  0.33333333,  0.33333333,  0.26666667,
                    0.        ,  0.26666667,  0.26666667,  0.13333333,  0.33333333),
                  c(0.46666667,  0.6       ,  0.33333333,  0.13333333,  0.13333333,
                    0.13333333,  0.13333333,  0.13333333,  0.53333333,  0.33333333,
                    0.26666667,  0.33333333,  1.        ,  0.73333333,  0.6       ,
                    0.        ,  0.13333333,  0.26666667,  0.26666667,  0.33333333),
                  c(0.26666667,  0.13333333,  0.26666667,  0.4       ,  0.2       ,
                    0.26666667,  0.33333333,  0.46666667,  0.4       ,  0.4       ,
                    0.4       ,  0.33333333,  0.73333333,  1.        ,  0.86666667,
                    0.13333333,  0.26666667,  0.33333333,  0.13333333,  0.53333333),
                  c(0.13333333,  0.13333333,  0.13333333,  0.26666667,  0.2       ,
                    0.4       ,  0.2       ,  0.46666667,  0.4       ,  0.26666667,
                    0.26666667,  0.26666667,  0.6       ,  0.86666667,  1.        ,
                    0.        ,  0.26666667,  0.2       ,  0.26666667,  0.4       ),
                  c(0.        ,  0.        ,  0.        ,  0.13333333,  0.        ,
                    0.        ,  0.        ,  0.        ,  0.26666667,  0.        ,
                    0.13333333,  0.        ,  0.        ,  0.13333333,  0.        ,
                    1.        ,  0.        ,  0.        ,  0.        ,  0.26666667),
                  c(0.26666667,  0.13333333,  0.13333333,  0.        ,  0.        ,
                    0.2       ,  0.        ,  0.        ,  0.06666667,  0.13333333,
                    0.        ,  0.26666667,  0.13333333,  0.26666667,  0.26666667,
                    0.        ,  1.        ,  0.2       ,  0.2       ,  0.        ),
                  c(0.4       ,  0.4       ,  0.53333333,  0.4       ,  0.26666667,
                    0.4       ,  0.33333333,  0.46666667,  0.46666667,  0.26666667,
                    0.4       ,  0.26666667,  0.26666667,  0.33333333,  0.2       ,
                    0.        ,  0.2       ,  1.        ,  0.46666667,  0.46666667),
                  c(0.13333333,  0.4       ,  0.66666667,  0.66666667,  0.53333333,
                    0.6       ,  0.46666667,  0.26666667,  0.6       ,  0.53333333,
                    0.4       ,  0.13333333,  0.26666667,  0.13333333,  0.26666667,
                    0.        ,  0.2       ,  0.46666667,  1.        ,  0.4       ),
                  c(0.6       ,  0.2       ,  0.2       ,  0.4       ,  0.13333333,
                    0.4       ,  0.26666667,  0.        ,  0.13333333,  0.33333333,
                    0.33333333,  0.33333333,  0.33333333,  0.53333333,  0.4       ,
                    0.26666667,  0.        ,  0.46666667,  0.4       ,  1.        ))
rownames(netmat3) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
colnames(netmat3) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
netmat3
class(netmat3)

isSymmetric(netmat3)  #TRUE

net3<- graph_from_adjacency_matrix( netmat3, mode = "undirected",weighted = TRUE, 
                                    diag = FALSE, add.colnames = NULL, add.rownames = NA) 
net3
class(net3)
laynet<-layout_in_circle(net3)
library(plotrix)
cs<-color.scale(seq(0,1,1/16), cs1=c(0,1),cs2=c(0,1),cs3=c(0,1), color.spec="rgb")
E(net3)$color<-cs[E(net3)$weight]
plot(net3, layout=laynet, vertex.frame.color="royalblue", vertex.color="royalblue",
     vertex.size=2.5, edge.width=(E(net3)$weight)*5, edge.color=E(net3)$color)
title(main=c("Normalized Netwrok from all previous matrices"))





#εξαετια cross correlation τ=1
netmat9 <- rbind(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0),
                 c(0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0),
                 c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0),
                 c(1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0),
                 c(0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                 c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0),
                 c(0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 c(1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                 c(1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0),
                 c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
rownames(netmat9) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
colnames(netmat9) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                       "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                       "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                       "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                       "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )

isSymmetric(netmat9)  #true

net9<- graph_from_adjacency_matrix( netmat9, mode = "undirected",weighted = NULL, 
                                    diag = FALSE, add.colnames = NULL, add.rownames = NA) 
graph.density(net9) #0.8421053
laynet<-layout_in_circle(net9)
plot(net9, layout=laynet, vertex.frame.color="royalblue", vertex.color="royalblue",
     vertex.size=2.5, edge.width=1.5, edge.arrow.size=0.5)
title(main=c("Network from pearson correlation with prewhitening"))




################################################################3
##normalized 6year windows from the first 3 2year windows
netmat12 <- rbind(c(1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                  c(0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0),
                  c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1),
                  c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1),
                  c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0),
                  c(0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0),
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1),
                  c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
                  c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1),
                  c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1),
                  c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1),
                  c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1),
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1),
                  c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0),
                  c(0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1))
rownames(netmat12) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                        "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                        "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                        "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                        "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
colnames(netmat12) <- rownames(netmat12)
  

isSymmetric(netmat12)  #false

net12<- graph_from_adjacency_matrix( netmat12, mode = "directed",weighted = NULL, 
                                     diag = FALSE, add.colnames = NULL, add.rownames = NA) 
graph.density(net12) #0.7552632
laynet<-layout_in_circle(net12)
plot(net12, layout=laynet, vertex.frame.color="royalblue", vertex.color="royalblue",
     vertex.size=2.5, edge.width=0.5, edge.arrow.size=0.5)
title(main=c("Directed network from MutInf for the first 6 years"))


###########################################################3
#granger causality
netmat12 <- rbind(   c(0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0),
                     c(0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0),
                     c(0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0),
                     c(0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0))
colnames(netmat12) <- c("Cocoa","CoffeeArabica", "CoffeeRobusta", "TeaAvg-3auctions", 
                        "TeaColombo", "TeaKolkata", "TeaMombasa", "Coconutoil", 
                        "Groundnutoil", "Soybeans", "Soybeanoil", "Soybeanmeal", 
                        "WheatCanadian", "WheatUS-SRW",  "WheatUS-HRW", "BananaUS", 
                        "Orange", "SugarEC" ,"SugarUS", "SugarWorld" )
rownames(netmat12) <- rev(colnames(netmat12))


isSymmetric(netmat12)  #false

net12<- graph_from_adjacency_matrix( netmat12, mode = "directed",weighted = NULL, 
                                     diag = FALSE, add.colnames = NULL, add.rownames = NA) 
graph.density(net12) #0.07894
laynet<-layout_in_circle(net12)
plot(net12, layout=laynet, vertex.frame.color="royalblue", vertex.color="royalblue",
     vertex.size=2.5, edge.width=0.4, edge.arrow.size=0.5)
title(main=c('Normalized network from CGCIM(1)"))
