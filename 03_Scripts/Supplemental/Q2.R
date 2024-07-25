# Figure S8 and S9

q1 <- read.csv("02_Data/results-noise.csv", sep = ";", dec = "." ) #read.csv("data/Q1.csv", sep = ";", dec = "," )
q <- q1[q1$treatment == "combined",]
q$value.cv <- as.numeric(q$value.cv)
q$value.prod <- as.numeric(q$value.prod)
head(q)
str(q)
  

colores1 <- c( "#e41a1c" ,"#377eb8", "#4daf4a", "#984ea3", "#ff7f00")


  
  #plot the main trends for fish and for enrichment

  
  for (i in 1:2){
    
    if(i ==1) { colores <- colores1[1:3] } else { colores <- colores1}
    if(i ==1) { col.reps <- rep(colores, each = 50) } else { col.reps <- rep(colores, each = 150) }
    
    if(i == 1) { q$cat <- q$nutrient_input } else { q$cat <- q$pop_n }
    if(i == 1) { q$cont <- q$abiotic } else { q$cont <- q$biotic }
    #for plots
    if( i == 1) { plot.text  <-  c("low", "med", "high") } else {
      plot.text <- c("8", "16", "32", "64", "128") }
    
    q.a <- q[q$part == "Aboveground",]
    q.t <- q[q$part == "Total",]
    
    #ag sup parts
    q.aa <- q.a[q.a$measure == "alpha",]
    q.ag <- q.a[q.a$measure == "gamma",]
    q.ab <- q.a[q.a$measure == "beta",]
    #tl sub parts
    q.ta <- q.t[q.t$measure == "alpha",]
    q.tg <- q.t[q.t$measure == "gamma",]
    q.tb <- q.t[q.t$measure == "beta",]
    
    #ag sums
    if(i ==1) { alpha.t <- tapply(q.ta$value.cv, list(q.ta$cat, q.ta$treatment), mean)[c(2,3,1),]
    gamma.t <- tapply(q.tg$value.cv, list(q.tg$cat, q.tg$treatment), mean)[c(2,3,1),]
    beta.t <- tapply(q.tb$value.cv, list(q.tb$cat, q.tb$treatment), mean)[c(2,3,1),] 
    
    alpha.a <- tapply(q.aa$value.cv, list(q.aa$cat, q.aa$treatment), mean)[c(2,3,1),]
    gamma.a <- tapply(q.ag$value.cv, list(q.ag$cat, q.ag$treatment), mean)[c(2,3,1),]
    beta.a <- tapply(q.ab$value.cv, list(q.ab$cat, q.ab$treatment), mean)[c(2,3,1),] } else {
      
      alpha.t <- tapply(q.ta$value.cv, list(q.ta$cat, q.ta$treatment), mean)[1:5,]
      gamma.t <- tapply(q.tg$value.cv, list(q.tg$cat, q.tg$treatment), mean)[1:5,]
      beta.t <- tapply(q.tb$value.cv, list(q.tb$cat, q.tb$treatment), mean)[1:5,] 
      
      alpha.a <- tapply(q.aa$value.cv, list(q.aa$cat, q.aa$treatment), mean)[1:5,]
      gamma.a <- tapply(q.ag$value.cv, list(q.ag$cat, q.ag$treatment), mean)[1:5,]
      beta.a <- tapply(q.ab$value.cv, list(q.ab$cat, q.ab$treatment), mean)[1:5,] }
    
    enrich.colors <- data.frame(col.reps = colores1[1:3],
                                nutrient_input = c("low", "medium", "high"))
    fish.colors <- data.frame(col.reps = colores1,
                              pop_n = c(8,16,32,64,128))
    #for TA
    if(i ==1) { q.ta <- merge(q.ta, enrich.colors)  } else 
    { q.ta <- merge(q.ta, fish.colors) }
    
    if(i ==1) { q.tg <- merge(q.tg, enrich.colors)  } else 
    { q.tg <- merge(q.tg, fish.colors) }
    
    if(i ==1) { q.tb <- merge(q.tb, enrich.colors)  } else 
    { q.tb <- merge(q.tb, fish.colors) }
    
    #for AA
    if(i ==1) { q.aa <- merge(q.aa, enrich.colors)  } else 
    { q.aa <- merge(q.aa, fish.colors) }
    
    if(i ==1) { q.ag <- merge(q.ag, enrich.colors)  } else 
    { q.ag <- merge(q.ag, fish.colors) }
    
    if(i ==1) { q.ab <- merge(q.ab, enrich.colors)  } else 
    { q.ab <- merge(q.ab, fish.colors) }
    
    #plots
    x11(width = 11, height = 6)
    par(mfrow = c(2,5))
    
    barplot( cbind(alpha.t, gamma.t), beside = T, 
             col = colores,
             main = ifelse(i == 1, "Enrichment Level TLPP", "Fish Pop TLPP"),
             names.arg = c("CValpha", "CVgamma"))
    legend(ifelse(i == 1, "topright","topright"), 
           plot.text, 
           text.col = colores,
           bty = "n")
    
    barplot( cbind(beta.t), beside = T, 
             col = colores,
             main = ifelse(i == 1, "Enrichment Level TLPP", "Fish Pop TLPP"),
             names.arg = c("CVbeta (PE)"))
    
    plot(q.ta$cont, q.ta$value.cv, col = q.ta$col.reps, main = "TLPP CValpha",
         xlab = ifelse(i == 1, "Enrichment Variation", "Fish Connectivity"),
         ylab = "CValpha")
    plot(q.tg$cont, q.tg$value.cv, col = q.tg$col.reps, main = "TLPP CVgamma",
         xlab = ifelse(i == 1, "Enrichment Variation", "Fish Connectivity"),
         ylab = "CVgamma")
    plot(q.tb$cont, q.tb$value.cv, col = q.tb$col.reps, main = 'TLPP PE',
         xlab = ifelse(i == 1, "Enrichment Variation", "Fish Connectivity"),
         ylab = "CVbeta (PE)")
    
    #Above ground
    barplot( cbind(alpha.a, gamma.a), beside = T, 
             col = colores,
             main = ifelse(i == 1, "Enrichment Level AGPP", "Fish Pop AGPP"),
             names.arg = c("CValpha", "CVgamma"))
    legend("topright", 
           plot.text, 
           text.col = colores,
           bty = "n")
    barplot( cbind( beta.a), beside = T, 
             col = colores,
             main = ifelse(i == 1, "Enrichment Level AGPP", "Fish Pop AGPP"),
             names.arg = c("CVbeta (PE)"))
    
    
    plot(q.aa$cont, q.aa$value.cv, col = q.aa$col.reps, main = "AGPP CValpha",
         xlab = ifelse(i == 1, "Enrichment Variation", "Fish Connectivity"),
         ylab = "CValpha")
    plot(q.ag$cont, q.ag$value.cv, col = q.ag$col.reps, main = "AGPP CVgamma",
         xlab = ifelse(i == 1, "Enrichment Variation", "Fish Connectivity"),
         ylab = "CVgamma")
    plot(q.ab$cont, q.ab$value.cv, col = q.ab$col.reps, main = 'AGPP PE',
         xlab = ifelse(i == 1, "Enrichment Variation", "Fish Connectivity"),
         ylab = "CVbeta (PE)")
    
    
    
  }
  
  
