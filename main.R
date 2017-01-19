# Tim Jak & Koen Veenenbos
# Team: JakVeenenbos
# 19 January 2017

# Libraries used, please download the RColorBrewer
library(raster)
library(sp)
library(randomForest)
library(hydroGOF)

# Loading data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")
load("data/trainingPoly.rda")


#Remove the water, clouds and cloud shadow pixel
GewataB1[GewataB1 > 750] = NA
GewataB2[GewataB2 > 1000] = NA
GewataB3[GewataB3 > 1250] = NA
GewataB4[GewataB4 < 1200 & GewataB4 > 4000] = NA
GewataB5[GewataB5 > 3250] = NA
GewataB7[GewataB7 > 2500] = NA
vcfGewata[vcfGewata > 100] = NA

# Brick all the data together
GewataAll = brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(GewataAll) = c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")
GewataAllDF = as.data.frame(getValues(GewataAll))

# Gives a summuary of GewataAllDF 
summary(GewataAllDF)

# Create plot set-up
par(mfrow=c(3,2))

# Plot the VCF versus the different bands 
# You could also make a function of this plotting method.
plot(VCF ~ band1, data = GewataAllDF, pch = ".", col = "chocolate1", xlim = c(0, 750), ylim = c(0, 100))
plot(VCF ~ band2, data = GewataAllDF, pch = ".", col = "chocolate1", xlim = c(0, 1000), ylim = c(0, 100))
plot(VCF ~ band3, data = GewataAllDF, pch = ".", col = "chocolate1", xlim = c(0, 1250), ylim = c(0, 100))
plot(VCF ~ band4, data = GewataAllDF, pch = ".", col = "chocolate1", xlim = c(0, 4000), ylim = c(0, 100))
plot(VCF ~ band5, data = GewataAllDF, pch = ".", col = "chocolate1", xlim = c(0, 3250), ylim = c(0, 100))
plot(VCF ~ band7, data = GewataAllDF, pch = ".", col = "chocolate1", xlim = c(0, 2500), ylim = c(0, 100))

# Creating the model
VCF_bands_lm = lm(VCF ~ band1 + band2 + band3 + band4 + band5 + band7, data=GewataAllDF)

# Squared value of the model answer
squared_lm = summary(VCF_bands_lm)$r.squared 

# Predict statment
pred_lm = predict(GewataAll, model = VCF_bands_lm, na.rm = TRUE)

# prepare plot
par(mfrow=c(1,2))

# Plot the prediction and the actual tree cover
plot(pred_lm, main = "Predicted tree cover")
plot(GewataAll$VCF, main = "Actual VCF tree cover")

# Root Mean Sqaure Error
RMSE = rmse(vcfGewata@data@values, pred_lm@data@values, na.rm=T)
#RMSE = sqrt(mean((vcfGewata@data@values - pred_lm@data@values)^2, na.rm=TRUE))

# Reproject the training polygons
repr_poly = spTransform(trainingPoly, proj4string(vcfGewata))

# Calculate RMSE per pixel
RMSE_pixel = sqrt(mean((vcfGewata - pred_lm)^2, na.rm=TRUE))

# Create a raster for the different classes 
RMSE_rast = rasterize(trainingPoly, RMSE_pixel, field = "Class")

# Create a zonal 
RMSE_zonal = zonal(RMSE_pixel, RMSE_rast, fun = "mean")

# Create a data frame for the rmse zonal
RMSE_zonal_df = as.data.frame(RMSE_zonal)

# Diffine the different classes
classes_RMSE = c("cropland", "forest", "wetland")

# Calculate rmse per class
RMSE_per_class = cbind(RMSE_zonal_df, classes_RMSE)

# Change names of the columns
names(RMSE_per_class) = c("Zone", "Mean RMSE", "Class")

