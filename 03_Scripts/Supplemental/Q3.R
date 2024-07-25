# Figure S13

q1 <- read.csv("02_Data/results-noise.csv", sep = ";", dec = "." )
q1$value.cv <- as.numeric(q1$value.cv)
q1$value.prod <- as.numeric(q1$value.prod)
head(q1)
str(q1)


q1.ag <- q1[q1$part == "Aboveground",]
q1.ag.a <- q1.ag[q1.ag$measure == "alpha",]
q1.ag.g <- q1.ag[q1.ag$measure == "gamma",]
q1.ag.b <- q1.ag[q1.ag$measure == "beta",]

q1.tl <- q1[q1$part == "Total",]
q1.tl.a <- q1.tl[q1.tl$measure == "alpha",]
q1.tl.g <- q1.tl[q1.tl$measure == "gamma",]
q1.tl.b <- q1.tl[q1.tl$measure == "beta",]


#hist(tq$value.prod)
#mean(tq$value.prod)
#median(tq$value.prod)
#hist(tq$value.cv)

#using only prod values greater than 6

x11(width = 12, height = 4)
par(mfrow = c(2,5), mar = c(3, 3.5, 2, 1.5),
    mgp = c(2,1,0))

for (i in 1:2) {
  jojo <- ifelse(i == 1, "Aboveground", "Total")
  no <- q1[q1$part == jojo,]
  
  now <- no[no$value.prod > mean(no$value.prod, na.rm = TRUE),]
  
  p.b <- now[now$measure == "beta",]
  p.a <- now[now$measure == "alpha",]
  p.g <- now[now$measure == "gamma",]
  
  # if(i == 1) { b.full <- q1.ag.b$value.cv} else {b.full <-  q1.tl.b$value.cv}
  # hist( b.full, main = paste("CV.a ", jojo ))
  # hist(p.b$value.cv, add = T, col = "gray")
  # 
  # if(i == 1) { a.full <- q1.ag.a$value.cv} else {a.full <-  q1.tl.a$value.cv}
  # hist( a.full, main = paste("CV.a ", jojo ))
  # hist(p.a$value.cv, add = T, col = "purple")
  # 
  # if(i == 1) { g.full <- q1.ag.g$value.cv} else {g.full <-  q1.tl.g$value.cv}
  # hist(g.full , main = paste("CV.g ", jojo ))
  # hist(p.g$value.cv, breaks = 15, add = T, col = "green")
  # 
  #max.val <- max(c(p.a$value.cv, p.g$value.cv))
  plot(p.a$value.prod~p.a$value.cv,
       #main = "CV.prod vs. CV.a",
       xlab = "CValpha",
       ylab = ifelse(i == 1, "AGPP", "TLPP") )
  plot(p.g$value.prod~p.g$value.cv, 
       #main = "CV.prod vs. CV.g",
       xlab = "CVgamma",
       ylab = ifelse(i == 1, "AGPP", "TLPP") )
  plot(p.b$value.prod~p.b$value.cv,
       #main = "CV.prod vs. CV.b",
       xlab = "CVbeta (PE)",
       ylab = ifelse(i == 1, "AGPP", "TLPP") )
  
  barplot(tapply(p.b$value.prod,  p.b$nutrient_input , mean ), beside = T , 
          xlab = "Nutrient enrichment",
          ylab = ifelse(i == 1, "AGPP", "TLPP"))
  barplot(tapply(p.b$value.prod,  p.b$pop_n , mean ), beside = T ,
          xlab = "Fish pop size",
          ylab = ifelse(i == 1, "AGPP", "TLPP"))
  #these show nothing
  #plot(p.b$value.prod~p.b$abiotic, main = "Enrichment var")
  #plot(p.b$value.prod~p.b$biotic, main = "Connectivity")
  
  
}


