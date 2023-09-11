################################################################################
################################################################################

dir.create(file.path("Figures"))

library(spdep)
library(INLA)
## Load data 

## Load data and SpatialPolygonsDataFrame
carto <- st_read("../Data/Carto/carto_gb.shp")
carto <- carto[order(carto$Code), ]

# Leemos la matriz de adyacencia (grafo conexo)
adj_gb <- as.matrix(read.table("../Data/Carto/adj_gb.txt"))


load("./Results_Model6c.Rdata")
res <- shared_rw1_t4


n <- length(unique(res$.args$data$Code))
t <- length(unique(res$.args$data$ID_time))


################################################################################
###########      Figure 7: Posterior median of the area-specific     ###########
###########    shared spatial effect and posterior median and 95%    ###########
###########     credible interval of the health-outcome-specific     ###########
###########                    temporal component                    ###########
################################################################################

############# posterior median of the area-specific shared spatial effect

marg.temporal.inc <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time1[[z]]))
inla.est.inc <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.inc[[z]]))


marg.spatial.inc <- lapply(1:n, function(z) inla.tmarginal(function(i) exp(i), res$marginals.random$ID_area[[z]]))
inla.est.inc <- lapply(1:n, function(z) inla.zmarginal(marg.spatial.inc[[z]]))
spat.inc <- sapply(1:n, function(z) round(inla.est.inc[[z]]$quant0.5,3))

marg.temporal.mort <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time2[[z]]))
inla.est.mort <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.mort[[z]]))


marg.spatial.mort <- lapply(n+1:n, function(z) inla.tmarginal(function(i) exp(i), res$marginals.random$ID_area[[z]]))
inla.est.mort <- lapply(1:n, function(z) inla.zmarginal(marg.spatial.mort[[z]]))
spat.mort <- sapply(1:n, function(z) round(inla.est.mort[[z]]$quant0.5,3))


ID_area <- unique(res$.args$data$Code)
data.spatial <- data.frame(cbind(c(ID_area, ID_area),c(spat.inc, spat.mort),
                                 c(rep("Incidence",n),rep("Mortality",n))))
colnames(data.spatial) <- c("ID_area","mean", "ID")
data.spatial<- data.spatial[order(data.spatial$ID, data.spatial$ID_area),]

data.spatial$mean <- as.numeric(data.spatial$mean)


library(dplyr)
df2<-left_join(data.spatial, carto, by = c('ID_area' = 'Code'))

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")

break_points <- c(round(seq(round(min(df$mean)-0.01,2), 1, length.out = 5),2),
                  round(seq(1, round(max(df$mean)+0.01,2), length.out = 5),2)[-1])

cols <- c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61",
          "#f46d43", "#d73027")


library(ggplot2)
s <- ggplot(df) +
  geom_sf(aes(fill = mean))+
  theme_void()+
  labs(title = 'Leukaemia cancer',
       subtitle = "",
       x = NULL,
       y = NULL) +
  facet_wrap(~ID) +
  theme(
    plot.title = element_text(size = 33, face = "bold",
                              hjust=0.5,color="gray35"),
    strip.text.x = element_text(size=25,hjust=0.1,vjust=0, face = "bold",
                                color="gray35"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    legend.position = 'bottom',
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.title = element_text(size = 20, color="gray35"),
    legend.text = element_text(size = 15, color="gray35"),
    legend.key = element_rect()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
               guide = "bins",
               breaks = break_points)+
  guides(fill = guide_colourbar(direction = 'horizontal',  
                                title=' ',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=5,
                                barwidth = 50,
                                barheight = 1))


############# posterior median and 95% credible interval of the health-outcome-
############# specific temporal component

marg.temporal.inc <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time1[[z]]))
inla.est.inc <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.inc[[z]]))
alpha.I <- exp(res$summary.fixed$mean[1])*10^5

data.temporal.inc <- data.frame(cbind(1:9 , rep("I",t), 
                                      sapply(1:t, function(z) round(inla.est.inc[[z]]$mean,3)),
                                      sapply(1:t, function(z) round(inla.est.inc[[z]]$quant0.5,3)),
                                      sapply(1:t, function(z) round(inla.est.inc[[z]]$quant0.025,3)),
                                      sapply(1:t, function(z) round(inla.est.inc[[z]]$quant0.975,3))))


marg.temporal.mort <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time2[[z]]))
inla.est.mort <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.mort[[z]]))
alpha.M <- exp(res$summary.fixed$mean[2])*10^5

data.temporal.mort <- data.frame(cbind(1:9 , rep("M",t), 
                                       sapply(1:t, function(z) round(inla.est.mort[[z]]$mean,3)),
                                       sapply(1:t, function(z) round(inla.est.mort[[z]]$quant0.5,3)),
                                       sapply(1:t, function(z) round(inla.est.mort[[z]]$quant0.025,3)),
                                       sapply(1:t, function(z) round(inla.est.mort[[z]]$quant0.975,3))))

data.temporal <- rbind(data.temporal.inc,data.temporal.mort)
colnames(data.temporal) <- c("ID_time","type","mean","median","inf","sup")
data.temporal$mean <- as.numeric(data.temporal$mean)
data.temporal$median <- as.numeric(data.temporal$median)
data.temporal$inf <- as.numeric(data.temporal$inf)
data.temporal$sup <- as.numeric(data.temporal$sup)


data.temporal$median <- c(alpha.I*data.temporal$median[1:t],alpha.M*data.temporal$median[(t+1):(2*t)])
data.temporal$inf <- c(alpha.I*data.temporal$inf[1:t],alpha.M*data.temporal$inf[(t+1):(2*t)])
data.temporal$sup <- c(alpha.I*data.temporal$sup[1:t],alpha.M*data.temporal$sup[(t+1):(2*t)])


library(ggplot2)

p <- ggplot(data=data.temporal, aes(x=ID_time, y=median, group=type)) + 
  geom_point(aes(shape=type, color=type), size = 3) +
  geom_line(aes(color=type), linewidth = 0.8) +
  scale_shape_manual(values=c(17,16), labels = c("Incidence","Mortality"))+
  scale_color_manual(values=c('#4575b4','#d73027'), labels = c("Incidence","Mortality")) +
  geom_ribbon(aes(ymin=inf,ymax=sup, colour=type),alpha = 0.3, 
              fill = c(rep("#abd9e9",t),rep("#fdae61",t)),
              col = "grey", linetype = "dotted") +
  scale_x_discrete(limits = c("1", "2","3", "4", "5", "6", "7", "8", "9"),
                   labels=c("1" = "2002-2003", "2" = " ",
                            "3" = "2006-2007", "4" = " ", "5" = "2010-2011",
                            "6" = " ", "7" = "2014-2015", "8" = " ",
                            "9" = "2018-2019")) +
  geom_hline(yintercept=alpha.I, linetype="dashed", color = '#4575b4', size=0.7) + 
  annotate("text", x=9, y=alpha.I - 0.5, size=7, 
           label=as.expression(expression(paste("exp(", alpha[I] ,")"))), color = '#4575b4') + 
  geom_hline(yintercept=alpha.M, linetype="dashed", color = '#d73027', size=0.7) + 
  annotate("text", x=9, y=alpha.M - 0.5, size=7, 
           label=as.expression(expression(paste("exp(", alpha[M] ,")"))), color = '#d73027') + 
  theme(
    plot.title = element_text(size = 25, face = "bold",
                              hjust=0.5,color="gray35"),
    plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank(),
    legend.text = element_text(size = 17, color="gray35"),
    legend.justification=c(-0.2,1.1), legend.position=c(0,1),
    axis.title = element_text(size = 20, color="gray35"),
    axis.text = element_text(size = 17, color="gray35")
  ) +
  labs(title="Leukaemia cancer") +
  ylab(as.expression(expression( paste("exp (", alpha[d] + gamma[td], ")") ))) +
  xlab("Period")


pdf("./Figures/Figure7.pdf", height=13, width=24)
library(ggpubr)
ga1 <- ggarrange(NULL, p, 
                 ncol = 1, nrow = 2)

ggarrange(s, ga1, 
          ncol = 2, nrow = 1)

dev.off()




################################################################################
###########   Figure 8: Posterior medians and 95% credible interval  ###########
###########          of the spatio-temporal effect in Cardiff,       ###########
###########   Lothian (Edinburgh), NHS Liverpool CCG and NHS Nort    ###########
###########                    East London CCG.                      ###########
################################################################################
marg.temporal.inc <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time1[[z]]))
inla.est.inc <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.inc[[z]]))
alpha.I <- exp(res$summary.fixed$mean[1])*10^5
marg.spattemp.inc <- lapply(1:n, function(l) lapply(1:t, function(z) inla.tmarginal(function(i) alpha.I*exp(i), 
                                                                                    res$marginals.random$ID_area1[[((z-1)*n)+l]])))
inla.est.inc <- lapply(1:n, function(l) lapply(1:t, function(z) inla.zmarginal(marg.spattemp.inc[[l]][[z]])))


marg.temporal.mort <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time2[[z]]))
inla.est.mort <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.mort[[z]]))
alpha.M <- exp(res$summary.fixed$mean[2])*10^5
marg.spattemp.mort <- lapply(1:n,function(l) lapply(1:t, function(z) inla.tmarginal(function(i) alpha.M*exp(i), 
                                                                                    res$marginals.random$ID_area1[[((z-1)*n)+l+n*t]])))
inla.est.mort <- lapply(1:n, function(l) lapply(1:t, function(z) inla.zmarginal(marg.spattemp.mort[[l]][[z]])))

spattemp.inc.median <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.inc[[l]][[z]]$quant0.5,3)))
spattemp.mort.median <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.mort[[l]][[z]]$quant0.5,3)))

ID_area <- unique(res$.args$data$Code)
data.spattemp <- data.frame(cbind(rep(ID_area,2*t) , rep(rep(1:t, each = n),2),
                                  c(rep("Incidence",n*t),rep("Mortality",n*t))))
data.spattemp<- data.spattemp[order(data.spattemp$X3, data.spattemp$X2, data.spattemp$X1),]

spattemp.inc <- NULL
spattemp.mort <- NULL
for (z in 1:t) {
  spattemp.inc <- c(spattemp.inc,sapply(1:n, function(l) spattemp.inc.median[[l]][z]))
  spattemp.mort <- c(spattemp.mort,sapply(1:n, function(l) spattemp.mort.median[[l]][z]))
}

data.spattemp <- cbind(data.spattemp, c(spattemp.inc,spattemp.mort))

colnames(data.spattemp) <- c("ID_area","ID_time","ID","mean")
data.spattemp$mean <- as.numeric(data.spattemp$mean)

spattemp.inc.inf <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.inc[[l]][[z]]$quant0.025,3)))
spattemp.mort.inf <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.mort[[l]][[z]]$quant0.025,3)))

spattemp.inc.sup <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.inc[[l]][[z]]$quant0.975,3)))
spattemp.mort.sup <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.mort[[l]][[z]]$quant0.975,3)))

spattemp1.inc <- NULL
spattemp2.inc <- NULL
spattemp1.mort <- NULL
spattemp2.mort <- NULL
for (z in 1:t) {
  spattemp1.inc <- c(spattemp1.inc,
                     sapply(1:n, function(l) spattemp.inc.inf[[l]][z]))
  spattemp2.inc <- c(spattemp2.inc,
                     sapply(1:n, function(l) spattemp.inc.sup[[l]][z]))
  
  spattemp1.mort <- c(spattemp1.mort,
                      sapply(1:n, function(l) spattemp.mort.inf[[l]][z]))
  spattemp2.mort <- c(spattemp2.mort,
                      sapply(1:n, function(l) spattemp.mort.sup[[l]][z]))
}

data.spattemp <- cbind(data.spattemp, 
                       c(spattemp1.inc,spattemp1.mort), 
                       c(spattemp2.inc,spattemp2.mort))
colnames(data.spattemp) <- c("ID_area","ID_time","ID","mean","inf","sup")
data.spattemp$mean <- as.numeric(data.spattemp$mean)
data.spattemp$inf <- as.numeric(data.spattemp$inf)
data.spattemp$sup <- as.numeric(data.spattemp$sup)


###Selected areas: 
selec.areas <- c("W06000015", "S08000024", "E38000101",
                 "E38000255")
names <- c("Cardiff", "Lothian (Edinburgh)", "NHS Liverpool CCG", "NHS North East London CCG")
pos <- which(data.spattemp$ID_area=="W06000015" |
               data.spattemp$ID_area=="S08000024" |
               data.spattemp$ID_area=="E38000101" |
               data.spattemp$ID_area=="E38000255")

inter <- data.spattemp[pos,]
min_y <- round(min(inter$inf),2)-0.01
max_y <- round(max(inter$sup),2)+0.01


pdf("./Figures/Figure8.pdf", height=6, width=10)
for (l in 1:length(selec.areas)) {
  p1 <- ggplot(data=inter[which(inter$ID_area==selec.areas[l]),], 
               aes(x=ID_time, y=mean, group = ID, color = ID)) + 
    geom_line(size = 1) +
    geom_point(size = 4, aes(shape = ID)) +
    scale_color_manual(values=c('#4575b4','#d73027')) +
    scale_shape_manual(values=c(17, 19)) +
    geom_ribbon(aes(ymin=inf,ymax=sup),alpha = 0.3,
                fill = c(rep("#abd9e9",t),rep("#fdae61",t)),
                col = "grey", linetype = "dotted") +
    ylim(min_y, max_y) +
    scale_x_discrete(limits = c("1", "2","3", "4", "5", "6", "7", "8", "9"),
                     labels=c("1" = "2002-2003", "2" = " ",
                              "3" = "2006-2007", "4" = " ", "5" = "2010-2011",
                              "6" = " ", "7" = "2014-2015", "8" = " ",
                              "9" = "2018-2019")) +
    geom_hline(yintercept=alpha.I, linetype="dashed", color = '#4575b4', size=0.7) + 
    annotate("text", x=9, y=alpha.I + 1, size=7, 
             label=as.expression(expression(paste("exp(", alpha[I] ,")"))), color = '#4575b4') + 
    geom_hline(yintercept=alpha.M, linetype="dashed", color = '#d73027', size=0.7) + 
    annotate("text", x=9, y=alpha.M + 1, size=7, 
             label=as.expression(expression(paste("exp(", alpha[M] ,")"))), color = '#d73027') + 
    theme(
      plot.title = element_text(size = 25, face = "bold",
                                hjust=0.5,color="gray35"),
      plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.title = element_blank(),
      legend.text = element_text(size = 17, color="gray35"),
      legend.justification=c(-0.2,1.1), legend.position=c(0,1),
      axis.title = element_text(size = 20, color="gray35"),
      axis.text = element_text(size = 17, color="gray35")
    ) +
    labs(title=paste(names[l])) +
    ylab(" ") +
    xlab("Period")
  
  
  print(p1)
}
dev.off()



################################################################################
###########        Figure C2: Posterior median of the spatially      ###########
###########         unstructured random effect for mortality         ###########
################################################################################

marg.temporal.mort <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time2[[z]]))
inla.est.mort <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.mort[[z]]))

marg.spatial.unst <- lapply(1:n, function(z) inla.tmarginal(function(i) exp(i), res$marginals.random$ID_unst[[z]]))
inla.est.unst <- lapply(1:n, function(z) inla.zmarginal(marg.spatial.unst[[z]]))
spat.unst <- sapply(1:n, function(z) round(inla.est.unst[[z]]$quant0.5,3))

ID_area <- unique(res$.args$data$Code)
data.spatial <- data.frame(cbind(ID_area , spat.unst))
data.spatial<- data.spatial[order(data.spatial$ID_area),]

colnames(data.spatial) <- c("ID_area","mean")
data.spatial$mean <- as.numeric(data.spatial$mean)

library(dplyr)
df2<-left_join(data.spatial, carto, by = c('ID_area' = 'Code'))

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")

break_points <- c(round(seq(round(min(df$mean)-0.01,2), 1, length.out = 5),2),
                  round(seq(1, round(max(df$mean)+0.01,2), length.out = 5),2)[-1])

cols <- c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61",
          "#f46d43", "#d73027")

library(ggplot2)
s <- ggplot(df) +
  geom_sf(aes(fill = mean))+
  theme_void()+
  labs(title = 'Leukaemia cancer',
       subtitle = "",
       x = NULL,
       y = NULL) +
  theme(
    plot.title = element_text(size = 33, face = "bold",
                              hjust=0.5,color="gray35"),
    strip.text.x = element_text(size=25,hjust=0.1,vjust=0, face = "bold",
                                color="gray35"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    plot.margin = margin(0.8, 1.5, 0.5, 1.5, "cm"),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    legend.position = 'bottom',
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.title = element_text(size = 20, color="gray35"),
    legend.text = element_text(size = 15, color="gray35"),
    legend.key = element_rect()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
               guide = "bins",
               breaks = break_points)+
  guides(fill = guide_colourbar(direction = 'horizontal',  
                                title=as.expression(expression( paste("exp (", u[i], ")") )),  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=5,
                                barwidth = 35,
                                barheight = 1))


pdf("./Figures/FigureC2.pdf", height=14, width=7.5)
print(s)
dev.off()



################################################################################
###########          Figure C3: Spatio-temporal distribution          ###########
################################################################################
marg.temporal.inc <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time1[[z]]))
inla.est.inc <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.inc[[z]]))
marg.spattemp.inc <- lapply(1:n, function(l) lapply(1:t, function(z) inla.tmarginal(function(i) exp(i), 
                                                                                    res$marginals.random$ID_area1[[((z-1)*n)+l]])))
inla.est.inc <- lapply(1:n, function(l) lapply(1:t, function(z) inla.zmarginal(marg.spattemp.inc[[l]][[z]])))


marg.temporal.mort <- lapply(1:t, function(z) inla.tmarginal(exp, res$marginals.random$ID_time2[[z]]))
inla.est.mort <- lapply(1:t, function(z) inla.zmarginal(marg.temporal.mort[[z]]))
marg.spattemp.mort <- lapply(1:n,function(l) lapply(1:t, function(z) inla.tmarginal(function(i) exp(i), 
                                                                                    res$marginals.random$ID_area1[[((z-1)*n)+l+n*t]])))
inla.est.mort <- lapply(1:n, function(l) lapply(1:t, function(z) inla.zmarginal(marg.spattemp.mort[[l]][[z]])))

spattemp.inc.median <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.inc[[l]][[z]]$quant0.5,3)))
spattemp.mort.median <- lapply(1:n, function(l) sapply(1:t, function(z) round(inla.est.mort[[l]][[z]]$quant0.5,3)))


ID_area <- unique(res$.args$data$Code)
data.spattemp <- data.frame(cbind(rep(ID_area,2*t) , rep(rep(1:t, each = n),2),
                                  c(rep("Incidence",n*t),rep("Mortality",n*t))))
data.spattemp<- data.spattemp[order(data.spattemp$X3, data.spattemp$X2, data.spattemp$X1),]

spattemp.inc <- NULL
spattemp.mort <- NULL
for (z in 1:t) {
  spattemp.inc <- c(spattemp.inc,sapply(1:n, function(l) spattemp.inc.median[[l]][z]))
  spattemp.mort <- c(spattemp.mort,sapply(1:n, function(l) spattemp.mort.median[[l]][z]))
}

data.spattemp <- cbind(data.spattemp, c(spattemp.inc,spattemp.mort))

colnames(data.spattemp) <- c("ID_area","ID_time","ID","mean")
data.spattemp$mean <- as.numeric(data.spattemp$mean)


###############
## spatial
##############
library(dplyr)
df2<-left_join(data.spattemp, carto, by = c('ID_area' = 'Code'))

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")

break_points <- c(round(seq(round(min(df$mean)-0.01,2), 1, length.out = 5),2),
                  round(seq(1, round(max(df$mean)+0.01,2), length.out = 5),2)[-1])

cols <- c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61",
          "#f46d43", "#d73027")

library(ggplot2)
s1 <- ggplot(df[which(df$ID=="Incidence"),]) +
  geom_sf(aes(fill = mean))+
  theme_void()+
  labs(title = 'Leukaemia cancer',
       subtitle = "Incidence",
       x = NULL,
       y = NULL) +
  facet_wrap(~ID_time, nrow = 2,labeller = as_labeller(c("1" = "2002-2003", "2" = "2004-2005 ",
                                                        "3" = "2006-2007", "4" = "2008-2009 ", 
                                                        "5" = "2010-2011","6" = "2012-2013 ", 
                                                        "7" = "2014-2015", "8" = "2016-2017 ",
                                                        "9" = "2018-2019"))) +
  
  theme(
    plot.title = element_text(size = 40, face = "bold",
                              hjust=0.5,color="gray35"),
    plot.subtitle = element_text(size = 35,hjust=0,color="gray35"),
    strip.text.x = element_text(size=30,hjust=0.1,vjust=0, face = "bold",
                                color="gray35"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    plot.margin = margin(0.8, 0.5, 3.7, 0.5, "cm"),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    legend.position = 'none',
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.title = element_text(size = 25, color="gray35"),
    legend.text = element_text(size = 20, color="gray35"),
    legend.key = element_rect()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
               guide = "bins",
               breaks = break_points)+
  guides(fill = guide_colourbar(direction = 'horizontal',  
                                title=" ", 
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=5,
                                barwidth = 0,
                                barheight = 1))

s2 <- ggplot(df[which(df$ID=="Mortality"),]) +
  geom_sf(aes(fill = mean))+
  theme_void()+
  labs(title = ' ',
       subtitle = "Mortality",
       x = NULL,
       y = NULL) +
  facet_wrap(~ID_time, nrow = 2, labeller = as_labeller(c("1" = "2002-2003", "2" = "2004-2005 ",
                                              "3" = "2006-2007", "4" = "2008-2009 ", 
                                              "5" = "2010-2011","6" = "2012-2013 ", 
                                              "7" = "2014-2015", "8" = "2016-2017 ",
                                              "9" = "2018-2019"))) +
  theme(
    plot.title = element_text(size = 40, face = "bold",
                              hjust=0.5,color="gray35"),
    plot.subtitle = element_text(size = 35,hjust=0,color="gray35"),
    strip.text.x = element_text(size=30,hjust=0.1,vjust=0, face = "bold",
                                color="gray35"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    legend.position = 'bottom',
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.title = element_text(size = 25, color="gray35"),
    legend.text = element_text(size = 20, color="gray35"),
    legend.key = element_rect()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
               guide = "bins",
               breaks = break_points,
               limits = c(round(min(df$mean),1), round(max(df$mean),1)),
               show.limits = FALSE)+
  guides(fill = guide_colourbar(direction = 'horizontal',  
                                title= ' ',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=5,
                                barwidth = 70,
                                barheight = 1))



pdf("./Figures/FigureC3.pdf", height=28, width=15)
library(ggpubr)
ggarrange(s1, s2,
          ncol = 1, nrow = 2)
dev.off()


################################################################################
###########    Figure C4: Posterior median of the evolution of the   ###########
########### geographical distribution of rates per 100,000 inhabitants##########
################################################################################
ID_area <- unique(res$.args$data$Code)
data.rates <- data.frame(cbind(rep(ID_area, t) , rep(1:t, each = n)))
data.rates<- data.rates[order(data.rates$X2, data.rates$X1),]
data.rates <- rbind(data.rates, data.rates)
data.rates <- cbind(rep(c("Incidence","Mortality"), each = n*t), data.rates,
                    (res$summary.fitted.values$`0.5quant`)*10^5)
colnames(data.rates) <- c("ID","ID_area","ID_time", "final.rates")

library(dplyr)
df2<-left_join(data.rates, carto, by = c('ID_area' = 'Code'))

write_sf(df2, "./prueba.shp")
df <-read_sf("./prueba.shp")


break_points <- c(round(seq(round(min(df$fnl_rts),0), 
                            round(max(df$fnl_rts),0), length.out = 14),0))
cols <- c("#fcffa4", "#f2e661", "#fac62d", "#fca50a", "#f8850f", "#ed6925",
          "#dd513a", "#c73e4c", "#ae305c", "#932667", "#781c6d", "#5d126e",
          "#420a68", "#240c4f", "#0c0826")

library(ggplot2)
s1 <- ggplot(df[which(df$ID=="Incidence"),]) +
  geom_sf(aes(fill = fnl_rts))+
  theme_void()+
  labs(title = 'Leukaemia cancer',
       subtitle = "Incidence",
       x = NULL,
       y = NULL) +
  facet_wrap(~ID_time, nrow = 2,labeller = as_labeller(c("1" = "2002-2003", "2" = "2004-2005 ",
                                                        "3" = "2006-2007", "4" = "2008-2009 ", 
                                                        "5" = "2010-2011","6" = "2012-2013 ", 
                                                        "7" = "2014-2015", "8" = "2016-2017 ",
                                                        "9" = "2018-2019"))) +
  theme(
    plot.title = element_text(size = 40, face = "bold",
                              hjust=0.5,color="gray35"),
    plot.subtitle = element_text(size = 35,hjust=0,color="gray35"),
    strip.text.x = element_text(size=30,hjust=0.1,vjust=0, face = "bold",
                                color="gray35"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    plot.margin = margin(0.8, 0.5, 3.7, 0.5, "cm"),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    legend.position = 'none',
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.title = element_text(size = 25, color="gray35"),
    legend.text = element_text(size = 20, color="gray35"),
    legend.key = element_rect()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
               guide = "bins",
               breaks = break_points)+
  guides(fill = guide_colourbar(direction = 'horizontal',  
                                title=" ",  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=5,
                                barwidth = 0,
                                barheight = 1))

s2 <- ggplot(df[which(df$ID=="Mortality"),]) +
  geom_sf(aes(fill = fnl_rts))+
  theme_void()+
  labs(title = ' ',
       subtitle = "Mortality",
       x = NULL,
       y = NULL) +
  facet_wrap(~ID_time, nrow = 2, labeller = as_labeller(c("1" = "2002-2003", "2" = "2004-2005 ",
                                                          "3" = "2006-2007", "4" = "2008-2009 ", 
                                                          "5" = "2010-2011","6" = "2012-2013 ", 
                                                          "7" = "2014-2015", "8" = "2016-2017 ",
                                                          "9" = "2018-2019"))) +
  theme(
    plot.title = element_text(size = 40, face = "bold",
                              hjust=0.5,color="gray35"),
    plot.subtitle = element_text(size = 35,hjust=0,color="gray35"),
    strip.text.x = element_text(size=30,hjust=0.1,vjust=0, face = "bold",
                                color="gray35"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    plot.margin = margin(0.8, 0.5, 0.5, 0.5, "cm"),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank(),
    legend.position = 'bottom',
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.title = element_text(size = 25, color="gray35"),
    legend.text = element_text(size = 20, color="gray35"),
    legend.key = element_rect()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
               guide = "bins",
               breaks = break_points,
               limits = c(round(min(df$fnl_rts),0), round(max(df$fnl_rts),0)),
               show.limits = FALSE)+
  guides(fill = guide_colourbar(direction = 'horizontal',  
                                title=as.expression(expression( paste(r[itd], " per 100,000 inhabitants") )),  ##rename default legend
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=5,
                                barwidth = 70,
                                barheight = 1))


pdf("./Figures/FigureC4.pdf", height=28, width=15)
library(ggpubr)
ggarrange(s1, s2,
          ncol = 1, nrow = 2)
dev.off()
