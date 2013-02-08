#load libraries
require(mgcv)
require(ggplot2) #for some very nice plotting
require(ggmap) #for including maps
require(RColorBrewer)

##########
#input data and look at it
gmpe.data <- read.table("Data/all_data.dat", header = TRUE)
head(gmpe.data)

##########
#some exploratory plots
p <- ggplot(gmpe.data, aes(factor(mag), PGA1))
p + geom_boxplot()

p <- ggplot(gmpe.data, aes(factor(siteClass), PGA1))
p + geom_boxplot()

p <- ggplot(gmpe.data, aes(factor(eq.type), PGA1))
p + geom_boxplot()

p <- ggplot(gmpe.data, aes(Repi, PGA1, alpha = .1))
p + geom_point() + theme(legend.position="none")

##########
#plot the stations on a map of Japan
stations <- unique(gmpe.data[,c("stnm", "stlon", "stlat", "siteClass")])
#uncomment the next line if you'd like to get the map yourself
#jp = get_map(location = "japan", zoom = 5, source = "osm", color = "bw")
#open street map's server is often unavailable, so I'd suggest saving this map
#uncomment the next line if you'd like to save the map you just downloaded
#save(jp, file = "jp")
load("jp")
locs <- ggmap(jp)
locs <- locs + 
  geom_point(data = stations, mapping = aes(x = stlon, y = stlat))
locs

##########
#load and plot the boundary stations
load("bndry.stns")
locs <- ggmap(jp)
locs <- locs + 
  geom_point(data = bndry.stns, mapping = aes(x = stlon, y = stlat))
locs

##########
#load the stations without the boundary stations and assign cross-validation sets
load("cv.stations")

#we have 332 non-boundary stations, so we will divide these up into 9 groups of 33 and 1 group of 35
set.number <- c(rep(1:9, each = 33), rep(10, 35))

#randomize the order of the stations
set.seed(30102)
cv.stations <- cv.stations[sample(1:nrow(cv.stations), size = nrow(cv.stations)), ]

#randomize the order of the cv set numbers
set.seed(1021)
set.number <- set.number[sample(1:length(set.number), size = length(set.number))]
cv.stations$set.number <- set.number

#most recent 5 earthquakes will be the testing set, so let's remove these from our data right now
testing.eqs <- tail(unique(gmpe.data$event.num),5)
testing.data <- gmpe.data[gmpe.data$event.num %in% testing.eqs, ]
training.data <- gmpe.data[!(gmpe.data$event.num %in% testing.eqs), ]

#change to character for later
training.data$stnm <- as.character(training.data$stnm)
testing.data$stnm <- as.character(testing.data$stnm)
cv.stations$stnm <- as.character(cv.stations$stnm)

##########
#simple cross-validation function
cv.fun <- function(cv.data, formula) {
  residuals <- c()
  for(i in 1:10) {
    test.stns <- cv.stations[cv.stations$set.number == i, "stnm"]
    train <- cv.data[!(cv.data$stnm %in% test.stns), ]
    test <- cv.data[cv.data$stnm %in% test.stns, ]
    model.fit <- gam(formula, data = train)
    predictions <- predict(model.fit, newdata = test)
    residuals <- c(residuals, test$PGA1 - predictions)
  }
  return(residuals)
}

##########
#our model formulas
#before we get too fancy, why don't we just fit a simple linear model, and then we can
#compare this to the fancier models
formula.lm = PGA1 ~ Repi + mag + depth + AVS30 + factor(siteClass) + factor(eq.type)

#now we can get a bit more complicated 
formula.gam = PGA1 ~ te(Repi, mag, by=eq.type, bs="tp", id=1) + 
  s(stlon, stlat, by = eq.type, bs="tp", id=1) + s(AVS30, bs = "cr") + 
  factor(siteClass) + factor(eq.type)

##########
#model fitting and checking
model.fit.lm <- gam(formula.lm, data = training.data)
model.fit.gam <- gam(formula.gam, data = training.data, method = "REML")

#summary of the linear model
summary(model.fit.lm)

#diagnostic plots of the linear model
gam.check(model.fit.lm)

#summary of the gam model
summary(model.fit.gam)

#diagnostic plots of the gam model
gam.check(model.fit.gam)

##########
#model evaluation
#cross-validation error
sqrt(mean(cv.fun(training.data, formula.lm)^2))
sqrt(mean(cv.fun(training.data, formula.gam)^2))

#testing error
predict.lm <- predict(model.fit.lm, newdata = testing.data)
predict.gam <- predict(model.fit.gam, newdata = testing.data)

sqrt(mean((testing.data$PGA1 - predict.lm)^2))
sqrt(mean((testing.data$PGA1 - predict.gam)^2))

##########
#compare with standard GMPE
sqrt(mean(testing.data$resids^2))

##########
#plot residuals
t1 <- rbind(training.data, training.data)
t1$model <- rep(c("gam", "linear"), each = nrow(training.data))
t1$residuals <- c(model.fit.gam$residuals, model.fit.lm$residuals)

brewer.div <- colorRampPalette(brewer.pal(9, "RdBu"), interpolate = "spline")
cols <- (brewer.div(256))
cutoffs <- seq(-1.7,1.7, length.out=257)
brks <- seq(-1.7,1.7,length.out=10)
labs <- round(brks,2)

p.spatial <- ggmap(jp)
p.spatial <- p.spatial + 
  geom_point(data = t1, mapping = aes(x = stlon, y = stlat, colour = residuals), size=4) +
  scale_colour_gradientn(colours=cols, 
                         name="Residual", values = cutoffs, limits = c(-1.7,1.7),
                         labels=labs, breaks=brks, rescaler =
                           function(x,...) x, oob=identity, guide="colourbar") +
  xlab("Longitude") + ylab("Latitude") + 
  labs(title = "Model residuals") +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size=16, vjust=.5), 
        axis.title.x = element_text(size=16, vjust=.25)) +
  facet_wrap(~model)
p.spatial