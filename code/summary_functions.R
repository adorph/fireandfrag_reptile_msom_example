#Check correlation between observed and mean estimated species richness  
checkCor <- function(data, model){
  ndf<-extract.par(model, "Nsite")
  ndf$C <- data$C
  ret <- cor(ndf$mean, ndf$C)
  return(ret)
}  

cumulativeDet <- function(data, model){
  require(ggplot2)
  require(dplyr)
  
  # This function calculates probability detected when present and associated 95% CRIs for 1 to 50 days 
  ddp50 <- function (x) {  
    day <- data.frame (day = seq(1, 50))        
    return (mutate (day, mean = 1-(1-x$est)^day, low = 1-(1-x$low)^day, upp = 1-(1-x$upp)^day))
  }
  
  spec <- data$species.list
  ndf <- data.frame()
  for(i in 1:length(spec)){
    det <- extract.par(model, "^lp\\[")[i,2] #Extract mean
    luci <- cbind(extract.par(model, "^lp\\[")[i,4], extract.par(model, "^lp\\[")[i,8]) #Extract 95% CRIs
    ddp <- data.frame(cbind (est = plogis (det), ci  = plogis (luci)))
    names (ddp) = c("est", "low", "upp")
    values <- ddp50(ddp)
    values$species <- spec[i]
    ndf <- rbind(ndf, values)
    rm(det, luci, ddp, values)
  }
  
  nddp <- ndf %>% filter(day %in% 5)
  mns <- format(round(nddp$mean, 2), nsmall=2)
  
  ndf <- ndf %>% mutate(species, lab=recode(species,
                                            "Eastern Three Lined Skink" = paste("ETS : mu(5) =", mns[1]),
                                            "Pale flecked Garden Skink" = paste("BUR : mu(5) =", mns[2]),
                                            "Robust Ctenotus" = paste("CBP : mu(5) =", mns[3]),
                                            "Shrubland Morethia Skink" = paste("DUA : mu(5) =", mns[4]),
                                            "Eastern Grey Kangarooo" = paste("EGK : mu(5) =", mns[5]),
                                            "Long nosed Bandicoot" = paste("LNB : mu(5) =", mns[6]), 
                                            "Long nosed Potoroo" = paste("LNP : mu(5) =", mns[7]), 
                                            "Red necked Wallaby" = paste("RNW : mu(5) =", mns[8]),
                                            "Short beaked Echidna" = paste("SBE : mu(5) =", mns[9]),
                                            "Southern Brown Bandicoot" = paste("SBB : mu(5) =", mns[10]),
                                            "Swamp Antechinus" = paste("SWA : mu(5) =", mns[11]),
                                            "Swamp Rat" = paste("SWR : mu(5) =", mns[12]),
                                            "Swamp Wallaby" = paste("SWW : mu(5) =", mns[13]),
                                            "White footed Dunnart" = paste("WFD : mu(5) =", mns[14]) ))
  
  #plot a graph
  aa <- ggplot(data=ndf, aes(x=day, y=mean, group=species))+
    geom_line(size=1)+
    geom_vline(xintercept=c(5), linetype="dotted")+  #change 20 to the deployment time for your cameras 
    geom_ribbon(aes(x=day, ymin=low, ymax=upp), color="NA", alpha=0.4, fill="grey")+
    labs(x="Days", y="Probability detected when present")+
    facet_wrap(.~ lab, scales="free_x", ncol=3, labeller=label_wrap_gen()) +
    theme_classic()+
    theme(axis.title=element_text(size=12))+
    theme(axis.text=element_text(size=10, colour="black"))+
    scale_x_continuous(limits=c(0,50), breaks=seq(0,50, 5), expand = c(0,0))+
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .2), expand = c(0,0))
  
  oddp <- ndf %>% filter(day %in% 1)
  sumdp <- cbind(oddp$species, dailyMean=oddp$mean, dailyLow=oddp$low, dailyHigh=oddp$upp, 
                 cumMean = nddp$mean, cumLow=nddp$low, cumHigh=nddp$upp)
  
  return(list(detPlot=aa, detDays=sumdp))
}

#See baseline estimates of species-specific occupancy and detection
occ.det.sum <- function(data, model){
  #See baseline estimates of species-specific occupancy and detection
  occ <- model$sims.list$lpsi
  det <- model$sims.list$lp
  nspec <- data$nspec
  #This includes occupancy and detection estimates for all observed species
  psi.occ <- plogis(occ[,1:nspec]) 
  p.det  <- plogis(det[,1:nspec]) 
  (occ.matrix <- cbind(apply(psi.occ,2,mean),apply(psi.occ,2,sd),  data$spec.name.list))
  (det.matrix <- cbind(apply(p.det,2,mean),apply(p.det,2,sd), data$spec.name.list))
  return(list(occ.matrix, det.matrix))
}

#Create a list of estimated true occurrence for species at a site, with X and Y coords for that site for visualisation
spp.site <- function(data, model) {
  df <- read.csv("C:/Users/adorph/Documents/3_Casterton_SpatialHeterogeneity/statistical_analysis/msom/site_coords.csv", header = T, sep = ",", dec = ".")
  df <- df[order(df$site),]
  df <- (df %>% select(site=site, xcoord, ycoord))[1:107,]
  spdf <- as.data.frame(model$mean$z)
  colnames(spdf) <- data$species.list
  spdf <- cbind(df, spdf)
  spsd <- as.data.frame(model$sd$z)
  colnames(spsd) <- data$species.list
  spsd <- cbind(df, spsd)
  sp50 <- as.data.frame(model$q50$z)
  colnames(sp50) <- data$species.list
  sp50 <- cbind(df, sp50)
  spdf <- spdf %>% pivot_longer(cols=4:ncol(spdf), names_to="Species", values_to = "mean")
  spsd <- spsd %>% pivot_longer(cols=4:ncol(spsd), names_to="Species", values_to = "sd")
  sp50 <- sp50 %>% pivot_longer(cols=4:ncol(sp50), names_to="Species", values_to = "q50")
  out <- merge(spdf, spsd, by=c("site", "xcoord", "ycoord", "Species"))
  out <- merge(out, sp50, by=c("site", "xcoord", "ycoord", "Species"))
  return(out)
}

#Extract rows from results based on parameter name
extract.par <- function(model, extract){
  ndf <- setDT(as.data.frame(model$summary), keep.rownames = "parameter")
  ndf <- filter(ndf, grepl(extract, parameter))
  return(ndf)
}

#Extract species richness with site Id
extract.nsite <- function(model){
  df <- read.csv("C:/Users/adorph/Documents/3_Casterton_SpatialHeterogeneity/statistical_analysis/msom/site_coords.csv",header = T, sep = ",", dec = ".")
  df <- df[order(df$site),]
  df <- (df %>% select(site=site, xcoord, ycoord))[1:107,]
  nsite <- extract.par(model, "Nsite")
  out <- cbind(df, nsite)
  out <- out[,-4]
  return(out)
}

#Plot observed number of species against estimated number of species
plotobs.est <- function(data, model){
  ndf<-extract.par(model, "Nsite")
  ndf$C <- data$C
  g <- ggplot(ndf) +
    geom_abline(aes(intercept=0, slope=1)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, x=C), width=0, color="grey") +
    geom_point(aes(x=C, y=mean)) +
    coord_cartesian(xlim=c(0,8), ylim=c(0,8)) +
    labs(x="Observed Number of Species", y = "Estimated Number of Species") +
    theme_bw() +
    theme(panel.grid=element_blank(),
          axis.text = element_text(size=9),
          axis.title = element_text(size=10))
  return(g) 
}

#Plot estimated (Nsite) against observed species richness
plotobs.sites <- function(data, model){
  ndf <- extract.par(model, "Nocc.fs")
  plot(data$obs.occ, ndf[, 2], xlab="Observed number of occupied sites", ylab = "Estimated version of quantity", ylim=c(0,max(ndf[,2])), frame=F, pch=16)
  abline(0,1)
  segments(data$obs.occ, ndf[,4], data$obs.occ, ndf[, 8], col="grey")
}

# Community distribution of average occupancy and detection probability
plotcomm.dist <- function(model){
  #Average species detection probability is 0.100 and average species occurrence is 0.35
  par(mfrow = c(1,2))       # Fig. 11-16
  psi.sample <- plogis(rnorm(10^6, mean = model$mean$mu.lpsi, sd = model$mean$sd.lpsi))
  p.sample <- plogis(rnorm(10^6, mean = model$mean$mu.lp, sd = model$mean$sd.lp))
  hist(psi.sample, freq = F, breaks = 50, col = "grey", xlab = "Species occupancy probability", ylab = "Density", main = "")
  abline(v=mean(psi.sample), col="red", lwd=2)
  hist(p.sample, freq = F, breaks = 50, col = "grey", xlab = "Species detection probability", ylab = "Density", main = "")
  abline(v=mean(p.sample), col="red", lwd=2)
}

#Plot occupancy vs detection estimate
plotocc.det <- function(data, model){
  occ.det <- occ.det.sum(data, model)
  occ <- setDT(as.data.frame(occ.det[[1]]), keep.rownames = "Species")
  colnames(occ) <- c("Species", "Mean.Occupancy", "SD.Occupancy", "id")
  det <- setDT(as.data.frame(occ.det[[2]]), keep.rownames = "Species")
  colnames(det) <- c("Species", "Mean.Detection", "SD.Detection", "id")
  occdet <- merge(occ, det, by=c("Species", "id"))
  
  g <- ggplot(occdet) +
    geom_errorbar(aes(ymin=(Mean.Detection-SD.Detection), ymax=(Mean.Detection+SD.Detection), x=Mean.Occupancy), width=0.0, color="grey") +
    geom_errorbarh(aes(xmin=Mean.Occupancy-SD.Occupancy, xmax=Mean.Occupancy+SD.Occupancy, y=Mean.Detection), height=0, color="grey") +
    geom_point(aes(x=Mean.Occupancy, y=Mean.Detection)) +
    theme_bw() +
    labs(y="Detection Estimate", x="Occupancy Estimate") +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    theme(panel.grid=element_blank(),
          axis.text=element_text(size=9),
          axis.title=element_text(size=10))
  
  return(g)
}

# Visualize covariate mean relationships for the average species
plotcovar.relation <- function(data, model){
  o.temp <- seq(13, 40,, 500)
  o.date <- seq(1, 146,, 500)
  
  temp.pred <- as.matrix((o.temp - mean(data$temp)) / sd(data$temp))
  date.pred <- as.matrix((o.date - mean(data$date)) / sd(data$date))
  
  # Predict detection for temp and rain
  # Put all predictions into a single array
  tmp <- model$sims.list                # grab MCMC samples
  nsamp <- length(tmp$mu.lp)                   # number of mcmc samples
  predC <- array(NA, dim = c(500, nsamp, 2)) # "C" for 'community mean'
  
  for(i in 1:nsamp){
    predC[,i,1] <-  plogis(tmp$mu.lp[i] + tmp$mu.betalp1[i] * temp.pred + tmp$mu.betalp2[i] * (temp.pred^2))
    predC[,i,2] <-  plogis(tmp$mu.lp[i] + tmp$mu.betalp3[i] * date.pred + tmp$mu.betalp4[i] * (date.pred^2))
  }
  
  # Get posterior means and 95% CRIs and plot
  pmC <- apply(predC, c(1,3), mean)
  criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975), na.rm=T))
  
  par(mfrow=c(1,2), mar=c(4,4,1,1))
  plot(o.temp, pmC[,1], ylim=c(0,1), col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, xlab = "Temp", ylab = "Community mean detection")
  matlines(o.temp, t(criC[,,1]), col = "grey", lty = 1)
  plot(o.date, pmC[,2], ylim=c(0,1), col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, xlab = "Date", ylab = "Community mean detection")
  matlines(o.date, t(criC[,,2]), col = "grey", lty = 1)
}

#Plot regression coefficients
plotreg.coef <- function(data, model, mcmc){
  all20 <- as.matrix(mcmc)
  str(all20)                    # look at the MCMC output
  pm <- apply(all20, 2, mean)    # Get posterior means and 95% CRIs
  pm1 <- setDT(as.data.frame(pm), keep.rownames = "parameter")
  cri <- apply(all20, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
  cri1 <- as.data.frame(cri)
  setDT(cri1, keep.rownames = "CRI")
  
  plotcoef <- function(data, mu.betalp, specno, mytitle, ylab){
    g <- ggplot(data=data) +
      geom_errorbarh(aes(xmin=`2.5%`, xmax=`97.5%`, y=id, color=sig1), size=1, height=0.01) +
      geom_point(mapping=aes(x=pm, y=id)) +
      geom_vline(aes(xintercept=0), color="black") +
      geom_vline(aes(xintercept=mu.betalp[,2]), color="red") +
      geom_vline(aes(xintercept=mu.betalp[,4]), linetype="dashed", color="red") +
      geom_vline(aes(xintercept=mu.betalp[,8]), linetype="dashed", color="red") +
      scale_color_manual(values=c("black", "blue")) + 
      scale_y_continuous(breaks=1:specno, labels=tdet$spec, trans = "reverse") +
      coord_cartesian(xlim=c(-2.5,2.5)) +
      labs(title=mytitle) +
      theme_bw() + 
      theme(legend.position = "none",
            panel.grid = element_blank(), 
            axis.title = element_blank(),
            plot.title = element_text(size=10),
            axis.text = element_text(size=9))
    if(ylab==F){
      g <- g + 
        theme(axis.text.y = element_blank())
    }
    return(g)
  }
  
  # Temp linear (Fig. 11 � 20 left)
  tlpm <- filter(pm1, grepl("betalp1", parameter))
  tlpm$id <- 1:data$nspec
  tlpm$spec <- data$species.list
  tlcri <- cri1 %>%
    select(CRI, contains("betalp1")) %>%
    pivot_longer(cols=-CRI, names_to="parameter", values_to="values") %>%
    pivot_wider(names_from="CRI", values_from = "values")
  tdet <- merge(tlpm, tlcri, by=c("parameter"))
  tdet$sig1 <- (tdet[,5] * tdet[,6]) > 0
  mu.betalp1 <- extract.par(model, "mu.betalp1")
  
  #Effects of temp (quadratic) on detection
  tqpm <- filter(pm1, grepl("betalp2", parameter))
  tqpm$id <- 1:data$nspec
  tqpm$spec <- data$species.list
  tqcri <- cri1 %>%
    select(CRI, contains("betalp2")) %>%
    pivot_longer(cols=-CRI, names_to="parameter", values_to="values") %>%
    pivot_wider(names_from="CRI", values_from = "values")
  tqdet <- merge(tqpm, tqcri, by=c("parameter"))
  tqdet$sig1 <- (tqdet[,5] * tqdet[,6]) > 0
  mu.betalp2 <- extract.par(model, "mu.betalp2")
  
  # Effects of date (linear) on detection
  dlpm <- filter(pm1, grepl("betalp3", parameter))
  dlpm$id <- 1:data$nspec
  dlpm$spec <- data$species.list
  dlcri <- cri1 %>%
    select(CRI, contains("betalp3")) %>%
    pivot_longer(cols=-CRI, names_to="parameter", values_to="values") %>%
    pivot_wider(names_from="CRI", values_from = "values")
  dldet <- merge(dlpm, dlcri, by=c("parameter"))
  dldet$sig1 <- (dldet[,5] * dldet[,6]) > 0
  mu.betalp3 <- extract.par(model, "mu.betalp3")
  
  #Effects of date (quadratic) on detection
  dqpm <- filter(pm1, grepl("betalp4", parameter))
  dqpm$id <- 1:data$nspec
  dqpm$spec <- data$species.list
  dqcri <- cri1 %>%
    select(CRI, contains("betalp4")) %>%
    pivot_longer(cols=-CRI, names_to="parameter", values_to="values") %>%
    pivot_wider(names_from="CRI", values_from = "values")
  dqdet <- merge(dqpm, dqcri, by=c("parameter"))
  dqdet$sig1 <- (dqdet[,5] * dqdet[,6]) > 0
  mu.betalp4 <- extract.par(model, "mu.betalp4")
  
  p1 <- plotcoef(tdet, mu.betalp1, data$nspec, "Temperature (linear)", ylab=T)
  p2 <- plotcoef(tqdet, mu.betalp2, data$nspec, "Temperature (quadratic)", ylab=F)
  p3 <- plotcoef(dldet, mu.betalp3, data$nspec, "Date (linear)", ylab=T)
  p4 <- plotcoef(dqdet, mu.betalp4, data$nspec, "Date (quadratic)", ylab=F)
  
  pg <- plot_grid(p3, p4, p1, p2, rel_widths = c(1.6,1))
  xlab <- textGrob("Parameter Estimates", gp=gpar(fontsize=10))
  ga <- grid.arrange(arrangeGrob(pg, bottom=xlab))
  return(ga)
}

# Predict detection for temperature and rainfall for each of the observed species
plotdet.covar <- function(data, model, mcmc){
  require(ggplot2)
  require(tidyr)
  
  o.temp <- seq(13, 40,, 500)
  o.date <- seq(1, 146,, 500)
  
  temp.pred <- as.matrix((o.temp - mean(data$temp)) / sd(data$temp))
  date.pred <- as.matrix((o.date - mean(data$date)) / sd(data$date))
  
  all20 <- as.matrix(mcmc)
  str(all20)                    # look at the MCMC output
  pm <- apply(all20, 2, mean)    # Get posterior means and 95% CRIs
  pm1 <- setDT(as.data.frame(pm), keep.rownames = "parameter")
  cri <- apply(all20, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
  cri1 <- as.data.frame(cri)
  cri2 <- cri1 %>% pivot_longer(cols=2:length(cri1), names_to = "parameter", values_to = "cri")
  # Effects of temp (linear) on detection
  par(mfrow = c(2,2), cex.lab = 0.7, cex.axis = 0.7)
  
  #lp
  lppm <- filter(pm1, grepl("^lp\\[", parameter))
  # Temp linear
  tlpm <- filter(pm1, grepl("betalp1", parameter))
  # Temp quadratic
  tqpm <- filter(pm1, grepl("betalp2", parameter))
  #Date linear
  dlpm <- filter(pm1, grepl("betalp3", parameter))
  # Date (quadratic) on detection
  dqpm <- filter(pm1, grepl("betalp4", parameter))
  
  predS <- array(NA, dim = c(500, data$nspec, 3))   # covariate value x species x response, "S" for 'species'
  p.coef <- cbind(lp=lppm[,2], betalp1 =tlpm[,2], betalp2 = tqpm[,2], betalp3 = dlpm[,2], betalp4 = dqpm[,2])
  
  for(i in 1:data$nspec){          # Loop over 16 observed species
    predS[,i,1] <- plogis(p.coef[i,1] + p.coef[i,2] * temp.pred + p.coef[i,3] * temp.pred^2)     # p ~ date
    predS[,i,2] <- plogis(p.coef[i,1] + p.coef[i,4] * date.pred + p.coef[i,5] * date.pred^2) # p ~ duration
  }
  
  # Plots for detection probability and temperature and rainfall
  dtemp <- as.data.frame(cbind(o.temp, predS[,,1]))
  cnames <- names(data$spec.name.list)
  colnames(dtemp) <- c("Values", cnames)
  dtemp$Variable <- c("Temperature")
  dtemp <- dtemp %>%
    select("Values", "Variable", names(data$spec.name.list)) %>%
    pivot_longer(cols = -c("Values", "Variable"), names_to = "Species", values_to = "Detection Probability")
  
  ddate <- as.data.frame(cbind(o.date, predS[,,2]))
  colnames(ddate) <- c("Values", cnames)
  ddate$Variable <- c("Date")
  ddate <- ddate %>%
    select("Values", "Variable", names(data$spec.name.list)) %>%
    pivot_longer(cols = -c("Values", "Variable"), names_to = "Species", values_to = "Detection Probability")
  
  ndat <- rbind(dtemp, ddate)
  
  ggplot(data=ndat, group=Variable) +
    geom_line(aes(x=Values, y=`Detection Probability`, color=Species), lwd=1) +
    facet_wrap(. ~ Variable, scales="free_x", ncol=2, strip.position = "bottom") +
    ylim(0,1) +
    theme_bw() +
    theme(strip.text = element_text(size = 9, margin = margin(0.0,0,0.0,0, "cm")),
          strip.background = element_rect(fill="white", color=NA),
          strip.placement = "outside",
          axis.text = element_text(size=9),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          legend.text = element_text(size=9),
          legend.title = element_text(size=9),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,0.05,0,0), "cm"),
          panel.spacing.x = unit(5, "mm"))
}