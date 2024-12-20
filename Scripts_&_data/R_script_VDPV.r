library(doMC)
library(fields)
library(gdistance)
library(geometry)
library(ks)
library(lubridate)
library(maptools)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(sampSurf)
library(seraphim)
library(sf)
library(sp)
library(spatstat)

savingPlots = TRUE
showingPlots = FALSE

if (!file.exists("Outbreak_results/One_outbreak_ID.csv")) # code used to generate a fictive dataset from RRW simulations:
	{
		data0 = read.csv("RRW2_simu.csv", head=T) # to load the simulated spatially-annotated tree file
		data0 = data0[which(!data0[,"node2"]%in%data0[,"node1"]),] # to only keep the tip nodes of the tree
		data1 = data0[,c("endLon","endLat","endYear")]; colnames(data1) = c("longitude","latitude","collection_date")
		data1 = data1[order(data1[,"collection_date"]),]; startYear = min(data1[,"collection_date"])
		days = matrix(nrow=dim(data1)[1], ncol=1); colnames(days) = "days"
		for (i in 1:dim(data1)[1]) days[i,1] = round((data1[i,"collection_date"]-startYear)*365.25)
		data1 = cbind(data1, days); maxDays = max(data1[,"days"])
		write.csv(data1, "Outbreak_results/One_outbreak_ID.csv", row.names=F, quote=F)
	}

datasets = list.files("Outbreak_results")
datasets = gsub("\\.csv","",datasets[grepl("\\.csv",datasets)])
datasets = datasets[which(datasets!="One_outbreak_ID")]
# datasets = datasets[which(datasets=="One_outbreak_ID")]
for (h in 1:length(datasets))
	{
		dir.create(file.path(paste0("Outbreak_results/",datasets[h],"_dir")), showWarnings=F)
	}

template = crop(raster("Template_rast.tif"), extent(-30,140,-28,50))
cellSurfaces = raster("WorldPop1km.tif")
cell_size = sqrt(cellSurfaces) # cell heights (or widths) in meters
popDensity_008 = crop(raster("Human_popD.tif"), extent(-30,140,-28,50))
	# NOTE: the "Human_popD.tif" file is too heavy for GitHub and is available on request
	# or can be retrieved from the WorldPop website (https://hub.worldpop.org/geodata/summary?id=29692)
popDensity_008_log = popDensity_008; popDensity_008_log[] = log(popDensity_008_log[]+1)
popDensity_08 = raster::aggregate(popDensity_008, 10, fun=sum)
popDensity_08_log = popDensity_08; popDensity_08_log[] = log(popDensity_08_log[]+1)
popDensity_5 = raster::aggregate(popDensity_008, 50, fun=sum)
popDensity_5_log = popDensity_5; popDensity_5_log[] = log(popDensity_5_log[]+1)
friction_008 = crop(raster("Friction_2015.tif"), extent(-30,140,-28,50))
	# NOTE: the "Friction_2015.tif" file is too heavy for GitHub and is available on request
	# or can be retrieved from the publication of Weiss et al. (2018, Nature)
friction_008[is.na(popDensity_008[])] = NA
friction_08 = raster::aggregate(friction_008, 10, fun=mean)
friction_5 = raster::aggregate(friction_008, 50, fun=mean)
polyExtent = as(extent(xmin(template), xmax(template), ymin(template), ymax(template)), "SpatialPolygons")  
international_polygons = shapefile("NaturalEarth_shp/Admin_0_country_borders.shp")
international_polygons = crop(international_polygons, extent(template))
international_borders = shapefile("NaturalEarth_shp/Admin_0_boundary_lines.shp")
international_borders = crop(international_borders, extent(template))

# 1. Subsetting the data based on kernel densties

for (h in 1:length(datasets))
	{
		if (!file.exists(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv")))
			{
				data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
				maxDays = max(data1[,"days"])
				template = raster("Template_rast.tif")
				template = crop(template, extent(min(data1[,1])-3,max(data1[,1])+3,min(data1[,2])-3,max(data1[,2])+3))
				rast = template; rast[] = NA
				gridSize = c(rast@ncols, rast@nrows)
				xyMin = c(rast@extent@xmin, rast@extent@ymin)
				xyMax = c(rast@extent@xmax, rast@extent@ymax)
				H = Hpi(data1[,c("longitude","latitude")])
				kde = kde(data1[,c("longitude","latitude")], H=H, compute.cont=T, gridsize=gridSize, xmin=xyMin, xmax=xyMax)
				r = raster(kde); c = rasterToContour(r, levels=kde$cont["5%"])
				contourLevel = contourLevels(kde, prob=0.05); polygons = list()
				contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
				for (i in 1:length(contourLines)) polygons[[i]] = Polygon(cbind(contourLines[[i]]$x,contourLines[[i]]$y))
				ps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps))
				spdf = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
				writeOGR(spdf, dsn="./Outbreak_results/", layer=datasets[h], driver="ESRI Shapefile")
				if (showingPlots == TRUE)
					{
						plot(data1[,c("longitude","latitude")], axes=F, ann=F)
						r = raster(kde); c = rasterToContour(r, levels=kde$cont["30%"]); lines(c, col="gray30")
						r = raster(kde); c = rasterToContour(r, levels=kde$cont["10%"]); lines(c, col="red")
						r = raster(kde); c = rasterToContour(r, levels=kde$cont["5%"]); lines(c, col="orange")
						r = raster(kde); c = rasterToContour(r, levels=kde$cont["1%"]); lines(c, col="green3")
					}
				buffer = data1[which(data1[,"days"]==0),]; data2 = data1[which(data1[,"days"]==0),]
				dir.create(file.path(paste0("./Outbreak_results/",layer=datasets[h],"_dir/KDE_95_contours/")), showWarnings=F)
				for (i in 1:maxDays)
					{
						if (showingPlots == FALSE) cat(paste0("day ",i,"\n"))
						indices1 = which(data1[,"days"]==i)
						if (length(indices1) > 0)
							{
								if (dim(data2)[1] >= 5)
									{
										indices2 = c(); col_pt = NA
										kde = kde(buffer[,c("longitude","latitude")], H=H, compute.cont=T, gridsize=gridSize, xmin=xyMin, xmax=xyMax)
										r = raster(kde); contour = rasterToContour(r, levels=kde$cont["5%"])
										threshold = contourLevels(kde, 0.05); r[r[]<threshold] = NA; r[!is.na(r[])] = i; crs(r) = crs(rast)
										writeRaster(r, paste0("./Outbreak_results/",layer=datasets[h],"_dir/KDE_95_contours/KDE_95_contours_day_",i-1,".tif"), overwrite=T)
										writeOGR(contour, dsn=paste0("./Outbreak_results/",layer=datasets[h],"_dir/KDE_95_contours/"),
												 layer=paste0("KDE_95_contours_day_",i-1), driver="ESRI Shapefile")
										if (showingPlots == TRUE) plot(contour, main=paste0("day ",i), cex.main=0.9, col.main="gray30")
										for (j in 1:length(indices1))
											{
												point_in_polygon = FALSE
												pt.x = data1[indices1[j],"longitude"]
												pt.y = data1[indices1[j],"latitude"]
												for (k in 1:length(contour@lines[[1]]@Lines))
													{
														pol.x = contour@lines[[1]]@Lines[[k]]@coords[,1]
														pol.y = contour@lines[[1]]@Lines[[k]]@coords[,2]
														if (point.in.polygon(pt.x, pt.y, pol.x, pol.y) != 0)
															{
																point_in_polygon = TRUE
															}
													}
												if (point_in_polygon == FALSE)
													{
														indices2 = c(indices2, indices1[j])
														col_pt = "green3"
													}	else	{
														col_pt = "red"
													}
												if (showingPlots == TRUE)
													{
														points(data1[indices1[j],c("longitude","latitude")], col=col_pt)
													}
											}
										if (length(indices2) > 0)
											{
												data2 = rbind(data2, data1[indices2,])
											}
									}	else	{
										data2 = rbind(data2, data1[indices1,])
									}
								buffer = rbind(buffer, data1[indices1,])
							}
					}
				data2 = unique(data2)
				index = which(data2[,"days"]==min(data2[,"days"]))
				buffer = data2[index,]
				lines = paste(data2[index,"longitude"],data2[index,"latitude"],sep="_")
				for (d in 1:max(data2[,"days"]))
					{
						indices1 = which(data2[,"days"]==d)
						if (length(indices1) > 0)
							{
								for (i in 1:length(indices1))
									{
										line = paste(data2[indices1[i],"longitude"],data2[indices1[i],"latitude"],sep="_")
										if (sum(lines==line) == 0)
											{
												lines = c(lines, line); buffer = rbind(buffer, data2[indices1[i],])
											}
									}
							}
					}
				data2 = buffer
				write.csv(data2, paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), quote=F, row.names=F)
				for (i in 1:maxDays)
					{
						if (file.exists(paste0("Outbreak_results/",datasets[h],"_dir/KDE_95_contours/KDE_95_contours_day_",i,".tif")))
							{
								kde = raster(paste0("Outbreak_results/",datasets[h],"_dir/KDE_95_contours/KDE_95_contours_day_",i,".tif"))
								rast = merge(rast, kde)
							}
					}
				writeRaster(rast, paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.tif"))
			}
	}

# 2. Estimating and plotting the wavefront velocities

for (h in 1:length(datasets))
	{
		template = crop(raster("Template_rast.tif"), extent(-30,140,-28,50))
		crs_AAEAC = "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
			# Africa Albers Equal Area Conic projection (AAEAC); https://spatialreference.org/ref/esri/africa-albers-equal-area-conic/
		crs_Behrmann = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
			# world Behrmann projection (https://gis.stackexchange.com/questions/138287/what-is-the-proj4-for-world-behrmann-54017)

		crs_metricSystem = crs_AAEAC; template_proj = raster("Template_Afri.tif") # to be uncommented when working on an "African" outbreak
		# crs_metricSystem = crs_Behrmann; template_proj = raster("Template_Asia.tif") # to be uncommented when working on an "Asian" outbreak

		# 2.1. Estimating the wavefront velocities

		data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
		data1_WGS84 = data.frame(lon=as.numeric(data1[,1]), lat=as.numeric(data1[,2]))
		coordinates(data1_WGS84) = c("lon", "lat"); proj4string(data1_WGS84) = CRS("+init=epsg:4326")
		data1_AAEAC = spTransform(data1_WGS84, crs_metricSystem); data1[,1:2] = data1_AAEAC@coords[,1:2]
		data2 = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), head=T)
		distances_WGS84 = matrix(nrow=dim(data2)[1], ncol=1)
		for (i in 1:dim(distances_WGS84)[1])
			{
				pt1 = cbind(data2[1,1],data2[1,2]); pt2 = cbind(data2[i,1],data2[i,2])
				distances_WGS84[i,1] = rdist.earth(pt1, pt2, miles=F)
			}
		data2_WGS84 = data.frame(lon=as.numeric(data2[,1]), lat=as.numeric(data2[,2]))
		coordinates(data2_WGS84) = c("lon", "lat"); proj4string(data2_WGS84) = CRS("+init=epsg:4326")
		data2_AAEAC = spTransform(data2_WGS84, crs_metricSystem); data2[,1:2] = data2_AAEAC@coords[,1:2]
		distances_AAEAC = matrix(nrow=dim(data2)[1], ncol=1)
		for (i in 1:dim(distances_AAEAC)[1])
			{
				distances_AAEAC[i,1] = sqrt(((data2[i,1]-data2[1,1])^2)+((data2[i,2]-data2[1,2])^2))
			}	# plot(distances_WGS84, distances_AAEAC/1000)
				
			# 2.1.1. Defining the analysis mask
	
		H = Hpi(data2[,c("longitude","latitude")])
		kde = kde(data2[,c("longitude","latitude")], H=H, compute.cont=T, gridsize=c(1000,1000))
		rast1 = raster(kde); contour = rasterToContour(rast1, levels=kde$cont["5%"])
		threshold = min(raster::extract(rast1, data2[,c("longitude","latitude")]))
		ps_list = list(); c = 0
		for (i in 1:length(contour@lines))
			{
				for (j in 1:length(contour@lines[[i]]@Lines))
					{
						p = Polygon(contour@lines[[i]]@Lines[[j]]@coords)
						c = c+1; ps_list[[c]] = Polygons(list(p),c)
					}
			}
		sps = SpatialPolygons(ps_list)
		rast2 = mask(rast1, sps)
		mask = raster::resample(rast2, template_proj)
		
				# 2.1.2. Interpolation of first time of invasion
		
		tps_model = Tps(x=data2[,c("longitude","latitude")], Y=data2[,"days"])
		tps = interpolate(template_proj, tps_model)
		tps_mask = crop(mask(tps, mask), extent(template_proj))
		
				# 2.1.3. Estimating the friction map
		
					# 2.1.3.1. Measuring the local slope/friction using a 3x3 moving windows filter
		
		raster_resolution = mean(c(res(tps_mask)[1],res(tps_mask)[2]))
		f = matrix(1/raster_resolution, nrow=3, ncol=3)
		f[c(1,3,7,9)] = 1/(sqrt(2)*raster_resolution); f[5] = 0
		fun = function(x, ...)
			{
				sum(abs(x-x[5])*f)/8
			}
		friction = focal(tps_mask, w=matrix(1,nrow=3,ncol=3), fun=fun, pad=T, padValue=NA, na.rm=F)
		
					# 2.1.3.2. Smoothing the resulting friction surface using an average 11x11 cell filter
		
		myAvSize = 11
		friction_sd10 = focal(friction, w=matrix(1/(myAvSize^2),nrow=myAvSize,ncol=myAvSize), pad=T, padValue=NA, na.rm=F)
		
				# 2.1.4. Estimating the spread rate
		
		spreadRate_sd10 = ((1/(friction_sd10))/1000)*7
		meanV = mean(spreadRate_sd10[],na.rm=T)
		cat("\tMean wavefront velocity for ",datasets[h]," = ",meanV," km/week\n",sep="")
						# Mean wavefront velocity for AFP_NIE_SOS_7 = 18.38376 km/week
						# Mean wavefront velocity for AFP_all_NIE_JIS = 22.97064 km/week
						# Mean wavefront velocity for One_outbreak_ID = 3.063116 km/week

		# 2.2. Plotting the wavefront velocities

		data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
		data2 = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), head=T)
	
			# 2.2.1. Defining the analysis mask
	
		H = Hpi(data2[,c("longitude","latitude")])
		kde = kde(data2[,c("longitude","latitude")], H=H, compute.cont=T, gridsize=c(1000,1000))
		rast1 = raster(kde); contour = rasterToContour(rast1, levels=kde$cont["5%"])
		threshold = min(raster::extract(rast1, data2[,c("longitude","latitude")]))
		ps_list = list(); c = 0
		for (i in 1:length(contour@lines))
			{
				for (j in 1:length(contour@lines[[i]]@Lines))
					{
						p = Polygon(contour@lines[[i]]@Lines[[j]]@coords)
						c = c+1; ps_list[[c]] = Polygons(list(p),c)
					}
			}
		sps = SpatialPolygons(ps_list)
		rast2 = mask(rast1, sps)
		mask = raster::resample(rast2, template)
		
				# 2.2.2. Interpolation of first time of invasion
		
		tps_model = Tps(x=data2[,c("longitude","latitude")], Y=data2[,"days"])
		tps = interpolate(template, tps_model)
		tps_mask = crop(mask(tps, mask), extent(template))
		
				# 2.2.3. Estimating the friction map
		
					# 2.2.3.1. Measuring the local slope/friction using a 3x3 moving windows filter
		
		raster_resolution = mean(c(res(tps_mask)[1],res(tps_mask)[2]))
		f = matrix(1/raster_resolution, nrow=3, ncol=3)
		f[c(1,3,7,9)] = 1/(sqrt(2)*raster_resolution); f[5] = 0
		fun = function(x, ...)
			{ 
				sum(abs(x-x[5])*f)/8
			}
		friction = focal(tps_mask, w=matrix(1,nrow=3,ncol=3), fun=fun, pad=T, padValue=NA, na.rm=F)
		
					# 2.2.3.2. Smoothing the resulting friction surface using an average 11x11 cell filter
		
		myAvSize = 11
		friction_sd10 = focal(friction, w=matrix(1/(myAvSize^2),nrow=myAvSize,ncol=myAvSize), pad=T, padValue=NA, na.rm=F)
				
				# 2.2.4. Plotting the resulting rasters
		
		if (savingPlots == TRUE)
			{
				successiveKDEs = raster(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.tif"))
				successiveKDEs[1] = max(data1[,"days"]); successiveKDEs[2] = 0
				tps_mask_cropped = tps_mask; tps_mask_cropped[is.na(template[])] = NA
				tps_mask_cropped[tps_mask_cropped[]<0] = 0
				friction_sd10_cropped = friction_sd10; friction_sd10_cropped[is.na(template[])] = NA
				spreadRate_sd10_cropped = spreadRate_sd10
				spreadRate_sd10_cropped[is.na(template[])] = NA
				r = spreadRate_sd10_cropped; r[(!is.na(r[]))&(r[]>1)] = 1
				spreadRate_sd10_truncated = r
				pdf(paste0("Outbreak_results/",datasets[h],".pdf"), width=10.0, height=4.6)
				par(mfrow=c(2,2), mar=c(0.0,0.0,0.0,0.4), oma=c(0.4,0.0,0.6,1.6), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				cols1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121])
				cols2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(181)[21:121])
				colsP = cols1[(((data1[,"days"]-min(data1[,"days"]))/(max(data1[,"days"])-min(data1[,"days"])))*100)+1]
				plot(polyExtent, ann=F, axes=F, col=NA, border=NA)
				plot(popDensity_008_log, col=colorRampPalette(brewer.pal(9,"Greys"))(17)[2:10], box=F, axes=F, legend=F, add=T)
				plot(international_polygons, add=T, lwd=0.2, col=NA, border="gray30", lty=1)
				for (i in dim(data1)[1]:1)
					{
						points(data1[i,c("longitude","latitude")], pch=16, cex=0.5, col=colsP[i])
						points(data1[i,c("longitude","latitude")], pch=1, cex=0.5, col="gray30", lwd=0.2)
					}
				mtext("1. First invasion times", side=3, line=-11.2, at=80, cex=0.7, font=1, col="gray30")
				mtext("(all, in days)", side=3, line=-12.0, at=80, cex=0.7, font=1, col="gray30")
				rect(xmin(template), ymin(template), xmax(template), ymax(template), xpd=T, lwd=0.2, border="gray30")
				legendRast = raster(as.matrix(seq(0,max(data1[,"days"]),1)))
				plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.945,0.955,0.040,0.960),
				     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
				     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.40,0)))
				colsP = cols1[(((data2[,"days"]-min(data2[,"days"]))/(max(data2[,"days"])-min(data2[,"days"])))*100)+1]
				plot(polyExtent, ann=F, axes=F, col=NA, border=NA)
				plot(popDensity_008_log, col=colorRampPalette(brewer.pal(9,"Greys"))(17)[2:10], box=F, axes=F, legend=F, add=T)
				plot(international_polygons, add=T, lwd=0.2, col=NA, border="gray30", lty=1)
				for (i in dim(data2)[1]:1)
					{
						points(data2[i,c("longitude","latitude")], pch=16, cex=0.5, col=colsP[i])
						points(data2[i,c("longitude","latitude")], pch=1, cex=0.5, col="gray30", lwd=0.2)
					}
				mtext("2. First invasion times", side=3, line=-11.2, at=80, cex=0.7, font=1, col="gray30")
				mtext("(filtered, in days)", side=3, line=-12.0, at=80, cex=0.7, font=1, col="gray30")
				rect(xmin(template), ymin(template), xmax(template), ymax(template), xpd=T, lwd=0.2, border="gray30")
				legendRast = raster(as.matrix(seq(0,max(data2[,"days"]),1)))
				plot(legendRast, legend.only=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.945,0.955,0.040,0.960),
				     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
				     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.40,0)))
				plot(polyExtent, ann=F, axes=F, col=NA, border=NA)
				plot(popDensity_008_log, col=colorRampPalette(brewer.pal(9,"Greys"))(17)[2:10], box=F, axes=F, legend=F, add=T)
				plot(tps_mask_cropped, col=cols1, colNA=NA, legend=F, add=T)
				plot(international_polygons, add=T, lwd=0.2, col=NA, border="gray30", lty=1)
				points(data2[,c("longitude","latitude")], pch=3, cex=0.4, lwd=0.4, col="gray30")
				mtext("3. First invasion times", side=3, line=-11.2, at=80, cex=0.7, font=1, col="gray30")
				mtext("(interpolated, in days)", side=3, line=-12.0, at=80, cex=0.7, font=1, col="gray30")
				rect(xmin(template), ymin(template), xmax(template), ymax(template), xpd=T, lwd=0.2, border="gray30")
				plot(tps_mask_cropped, legend.only=T, add=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.945,0.955,0.040,0.960),
				     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
				     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.40,0)))
				plot(polyExtent, ann=F, axes=F, col=NA, border=NA)
				plot(popDensity_008_log, col=colorRampPalette(brewer.pal(9,"Greys"))(17)[2:10], box=F, axes=F, legend=F, add=T)
				plot(friction_sd10_cropped, col=cols1, colNA=NA, legend=F, add=T)
				plot(international_polygons, add=T, lwd=0.2, col=NA, border="gray30", lty=1)
				points(data2[,c("longitude","latitude")], pch=3, cex=0.4, lwd=0.4, col="gray30")
				mtext("4. Friction map", side=3, line=-11.2, at=80, cex=0.7, font=1, col="gray30")
				rect(xmin(template), ymin(template), xmax(template), ymax(template), xpd=T, lwd=0.2, border="gray30")
				plot(friction_sd10_cropped, legend.only=T, add=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.945,0.955,0.040,0.960), 
				     alpha=1, legend.args=list(text="", cex=0.5, line=0.5, col="gray30"), axis.args=list(cex.axis=0.6, lwd=0,
				     lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.40,0)))
				dev.off()
			}
	}
	
# 3. Investigating the impact of environmental factors on the dispersal

source("circuitscapeFct.r")
environmentalExtraction = function(segments, rasters, resistances, pathModels, ID="obs")
	{
		if (!is.list(rasters)) nberOfColumns = length(pathModels)
		if (is.list(rasters)) nberOfColumns = length(pathModels)*length(rasters)
		extractions = matrix(nrow=dim(segments)[1], ncol=nberOfColumns)
		fromCoor = matrix(nrow=dim(segments)[1], ncol=2)
		firstCase = t(segments[which(segments[,"d1"]==0)[1],c("x1","y1")])
		for (i in 1:dim(fromCoor)[1]) fromCoor[i,] = firstCase
		toCoor = segments[,c("x2","y2")]
		for (i in 1:length(pathModels))
			{
				if (pathModels[i] == 1)
					{
						for (j in 1:dim(segments)[1])
							{
								point1 = segments[which(segments[,"d1"]==0)[1],c("x1","y1")]
								point2 = segments[j,c("x2","y2")]; names(point2) = c("x","y")
								points = rbind(point1,point2); linesList = list()
								linesList[[1]] = Lines(list(Line(points)),1)
								spatialLine = SpatialLines(linesList)
								for (k in 1:length(rasters))
									{
										values = raster::extract(rasters[[k]], spatialLine)[[1]]
										if (resistances[k] == FALSE) values = 1/values
										column = ((i-1)*length(pathModels))+k
										extractions[j,column] = sum(values, na.rm=T)
									}
							}
					}
				if (pathModels[i] == 2)
					{
						for (j in 1:length(rasters))
							{
								if (resistances[j] == FALSE)
									{
										trEnvVariable = transition(rasters[[j]], mean, directions=8)
									}	else	{
										trEnvVariable = transition(rasters[[j]], function(x) 1/mean(x), directions=8)
									}
								trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
								envDistances = costDistance(trEnvVariableCorr, fromCoor[1,], toCoor)
								column = ((i-1)*length(pathModels))+j; extractions[,column] = t(envDistances)
							}
					}
				if (pathModels[i] == 3)
					{
						wd1 = getwd(); wd2 = paste0(wd1,"/Circuitscape_tmp"); setwd(wd2); buffer = list()
						for (j in 1:length(rasters))
							{
								rasterName = paste0(names(rasters[[j]]))
								if (!file.exists(paste0(rasterName,".asc")))
									{
										writeRaster(rasters[[j]], paste0(rasterName,".asc"), overwrite=T, showWarnings=F)
									}
								circuitscapeFct(rasters[[j]], rasterName, resistances[j], resistances[j], 
												fourCells=F, fromCoor, toCoor, OS="Unix", rasterName, ID)
							}
						for (j in 1:length(rasters))
							{
								rasterName = paste0(names(rasters[[j]]))
								folder = paste(rasterName,"_",ID,sep="")
								if (file.exists(paste0(folder,"/raster_file_temp_resistances.txt")))
									{
										tab = read.table(paste0(folder,"/raster_file_temp_resistances.txt"), header=F)
										tab = tab[2:dim(tab)[1], 2:dim(tab)[2]]
										mat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=1)
										sameCoordinates = FALSE
										if (sum(fromCoor-toCoor) == 0) sameCoordinates = TRUE
										if (sameCoordinates == FALSE)
											{
												for (k in 1:length(fromCoor[,1])) mat[k] = tab[k,(k+length(fromCoor[,1]))]
											}	else	{
												mat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=dim(as.matrix(fromCoor))[1])
												mat = tab[1:dim(fromCoor)[1],1:dim(fromCoor)[1]]
											}
										column = ((i-1)*length(pathModels))+j; extractions[,column] = mat
										unlink(folder, recursive=T)
									}	else	{
										column = ((i-1)*length(pathModels))+j; extractions[,column] = NA
									}
							}
						setwd(wd1)
					}
			}
		colNames = c()
		for (i in 1:length(pathModels))
			{
				if (pathModels[i] == 1) pathModel = "SL"
				if (pathModels[i] == 2) pathModel = "LC"
				if (pathModels[i] == 3) pathModel = "CS"
				for (j in 1:length(rasters)) colNames = c(colNames, paste0(pathModel,"_",names(rasters[[j]])))
			}
		colnames(extractions) = colNames
		return(extractions)
	}
rasters = list(); resistances = c(); c = 0
popDensity = popDensity_08_log
# popDensity = popDensity_5_log
friction = friction_08
null_raster = popDensity
names(null_raster) = "Null_raster"
null_raster[!is.na(null_raster[])] = 1
c = c + 1; rasters[[c]] = null_raster
resistances = c(resistances, TRUE)
kS = c(10, 100, 1000, 10000, 100000)
for (k in kS)
	{
		c = c + 1
		if (!file.exists(paste0("Prepared_rasters/Population_density_log_C_k",k,".tif")))
			{
				r = popDensity; r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
				names(r) = paste0("Population_density_log_C_k",k)
				writeRaster(r, paste0("Prepared_rasters/Population_density_log_C_k",k,".tif"))
			}	else	{
				r = raster(paste0("Prepared_rasters/Population_density_log_C_k",k,".tif"))
			}
		names(r) = paste0("Population_density_log_C_k",k)
		rasters[[c]] = r; resistances = c(resistances, FALSE)
	}
for (k in kS)
	{
		c = c + 1
		if (!file.exists(paste0("Prepared_rasters/Friction_R_k",k,".tif")))
			{
				r = friction; r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
				names(r) = paste0("Friction_R_k",k)
				writeRaster(r, paste0("Prepared_rasters/Friction_R_k",k,".tif"))
			}	else	{
				r = raster(paste0("Prepared_rasters/Friction_R_k",k,".tif"))
			}
		names(r) = paste0("Friction_R_k",k)
		rasters[[c]] = r; resistances = c(resistances, TRUE)
	}
for (k in kS)
	{
		c = c + 1
		if (!file.exists(paste0("Prepared_rasters/International_borders_R_k",k,".tif")))
			{
				r = rasterize(international_borders, null_raster, fun=min)
				r[!is.na(r[])] = 1; r[(is.na(r[]))&(!is.na(null_raster[]))] = 0
				r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
				names(r) = paste0("International_borders_R_k",k)
				writeRaster(r, paste0("Prepared_rasters/International_borders_R_k",k,".tif"))
			}	else	{
				r = raster(paste0("Prepared_rasters/International_borders_R_k",k,".tif"))
			}
		names(r) = paste0("International_borders_R_k",k)
		rasters[[c]] = r; resistances = c(resistances, TRUE)
	}
for (k in kS)
	{
		c = c + 1
		if (!file.exists(paste0("Prepared_rasters/International_borders_C_k",k,".tif")))
			{
				r = rasterize(international_borders, null_raster, fun=min)
				r[!is.na(r[])] = 1; r[(is.na(r[]))&(!is.na(null_raster[]))] = 0
				r[!is.na(r[])] = 1+(k*(r[!is.na(r[])]))
				names(r) = paste0("International_borders_C_k",k)
				writeRaster(r, paste0("Prepared_rasters/International_borders_C_k",k,".tif"))
			}	else	{
				r = raster(paste0("Prepared_rasters/International_borders_C_k",k,".tif"))
			}
		names(r) = paste0("International_borders_C_k",k)
		rasters[[c]] = r; resistances = c(resistances, FALSE)
	}

	# 3.1. Rotating the observed segments to generate a null dispersal model

nberOfRandomisations = 100
for (h in 1:length(datasets))
	{
		data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
		data2 = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), head=T)
		segments_obs = c(); index0 = which(data2[,"days"]==0)
		for (i in 1:dim(data2)[1])
			{
				if (i != index0)
					{
						segment = cbind(data2[index0,"longitude"],data2[index0,"latitude"])
						segment = cbind(segment, data2[i,"longitude"], data2[i,"latitude"])
						segment = cbind(segment, data2[index0,"days"], data2[i,"days"])
						segments_obs = rbind(segments_obs, cbind(segment, data2[i,"days"]-data2[index0,"days"]))
					}
			}
		colnames(segments_obs) = c("x1","y1","x2","y2","d1","d2","dt")
		write.csv(segments_obs, paste0("Outbreak_results/",datasets[h],"_dir/Obs_segments.csv"), row.names=F, quote=F)	
		dir.create(file.path(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/")), showWarnings=F)
		for (i in 1:nberOfRandomisations)
			{
				if (!file.exists(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_segments_",i,".csv")))
					{
						rotations = segments_obs; # print(c(i))
						for (j in 1:dim(rotations)[1])
							{
								inLand = FALSE; counter = 0
								while (inLand == FALSE)
									{
										x1 = segments_obs[j,"x1"]; y1 = segments_obs[j,"y1"]
										x2 = segments_obs[j,"x2"]; y2 = segments_obs[j,"y2"]
										angle = (2*pi)*runif(1); s = sin(angle); c = cos(angle)
										x = x2-x1; y = y2-y1
										x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
										x_new = x_new+x1; y_new = y_new+y1
										if (!is.na(raster::extract(rasters[[1]],cbind(x_new,y_new))))
											{
												inLand = TRUE
												rotations[j,"x2"] = x_new; rotations[j,"y2"] = y_new
											}	else	{
												counter = counter+1
												if (counter == 100)
													{
														inLand = TRUE; print(c(i,j,counter))
													}
											}
									}
							}		
						write.csv(rotations, paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_segments_",i,".csv"), row.names=F, quote=F)
					}
			}
	}
	
	# 3.2. Analysing the impact of factors on the wavefront dispersal velocity

nberOfRandomisations = 100; pathModels = c(3); pathModel = "CS"; nberOfCores = 10; registerDoMC(cores=nberOfCores)
for (h in 1:length(datasets))
	{
		data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
		data2 = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), head=T)
		segments_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_segments.csv")); ID = "obs"
		extractions_obs = environmentalExtraction(segments_obs, rasters, resistances, pathModels, ID)
		write.csv(extractions_obs, paste0("Outbreak_results/",datasets[h],"_dir/Obs_extraction.csv"), row.names=F, quote=F)
		extractions_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_extraction.csv"), head=T)
		rasters_to_select = c(1)
		for (i in 2:length(rasters))
			{
				envVariableName = names(rasters[[i]])
				extractions_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_extraction.csv"), header=T)
				rasterNull = paste0(pathModel,"_Null_raster")
				rasterName = paste0(pathModel,"_",envVariableName)
				y = segments_obs[,"d2"]; x = extractions_obs[,rasterNull]
				LR = lm(as.formula(paste0("y ~ x")))
				R2_obs_nul = summary(LR)$r.squared
				y = segments_obs[,"d2"]; x = extractions_obs[,rasterName]
				LR = lm(as.formula(paste0("y ~ x")))
				R2_obs_env = summary(LR)$r.squared
				Q_obs = R2_obs_env-R2_obs_nul
				if (Q_obs > 0) rasters_to_select = c(rasters_to_select, i)
			}
		rasters_selected = rasters[rasters_to_select]
		resistances_selected = resistances[rasters_to_select]
		buffer = foreach(i = 1:nberOfRandomisations) %dopar% {
		# for (i in 1:nberOfRandomisations) {
				segments_ran = as.matrix(read.csv(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_segments_",i,".csv"), head=T))
				extractions_ran = environmentalExtraction(segments_ran, rasters_selected, resistances_selected, pathModels, ID=paste0("ran",i))
				write.csv(extractions_ran, paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_extraction_",i,".csv"), row.names=F, quote=F)
				i
			}
	}
Qs_obs_list1 = list(); Qs_ran_list1 = list()
for (h in 1:length(datasets))
	{
		Qs_obs_list2 = list(); Qs_ran_list2 = list(); n = 0
		segments_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_segments.csv")); ID = "obs"
		for (i in 2:length(rasters))
			{
				for (j in 1:length(pathModels))
					{
						if (pathModels[j] == 2) { pathModel = "LC" }
						if (pathModels[j] == 3) { pathModel = "CS" }
						envVariableName = names(rasters[[i]])
						extractions_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_extraction.csv"), header=T)
						rasterNull = paste0(pathModel,"_Null_raster")
						rasterName = paste0(pathModel,"_",envVariableName)
						y = segments_obs[,"d2"]; x = extractions_obs[,rasterNull]
						LR = lm(as.formula(paste0("y ~ x")))
						R2_obs_nul = summary(LR)$r.squared
						y = segments_obs[,"d2"]; x = extractions_obs[,rasterName]
						LR = lm(as.formula(paste0("y ~ x")))
						R2_obs_env = summary(LR)$r.squared
						Q_obs = R2_obs_env-R2_obs_nul
						Q_rans = rep(NA, nberOfRandomisations); c = 0
						if (Q_obs <= 0)
							{
								cat(datasets[[h]]," - ",rasterName,", Qobs = ",round(Q_obs,3),"\n",sep="")
							}	else	{
								for (k in 1:nberOfRandomisations)
									{
										segments_ran = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_segments_",k,".csv"), header=T)
										extractions_ran = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_extraction_",k,".csv"))
										y = segments_ran[,"d2"]; x = extractions_ran[,rasterNull]
										LR = lm(as.formula(paste0("y ~ x")))
										R2_ran_nul = summary(LR)$r.squared
										y = segments_ran[,"d2"]; x = extractions_ran[,rasterName]
										LR = lm(as.formula(paste0("y ~ x")))
										R2_ran_env = summary(LR)$r.squared
										Q_ran = R2_ran_env-R2_ran_nul; Q_rans[k] = Q_ran
										if (Q_obs < Q_ran) c = c+1
									}
								pValue = round(c/nberOfRandomisations,3)
								cat(datasets[[h]]," - ",rasterName,", Qobs = ",round(Q_obs,3),", p-value = ",pValue,"\n",sep="")
							}
						n = n+1; Qs_obs_list2[[n]] = Q_obs; Qs_ran_list2[[n]] = Q_rans
					}
			}
		saveRDS(Qs_obs_list2, paste0("Outbreak_results/",datasets[h],"_dir/Q_obs_values.rds"))
		saveRDS(Qs_ran_list2, paste0("Outbreak_results/",datasets[h],"_dir/Q_ran_values.rds"))
		Qs_obs_list1[[h]] = Qs_obs_list2; Qs_ran_list1[[h]] = Qs_ran_list2
	}

	# 3.3. Analysis of the impact of barriers on the dispersal frequency

countingCrossingBarrierEvents = function(segments, shapefiles)
	{
		I = which(segments[,"d1"]==0)[1]
		crossingBarrierEvents = matrix(nrow=dim(segments)[1], ncol=length(shapefiles))
		if (showingPlots == TRUE)
			{
				pdf(paste0("TEMP.pdf")); plot(shapefiles[[1]], lwd=0.5, col="gray30")
			}
		for (i in 1:dim(segments)[1])
			{
				point1 = segments[I,c("x1","y1")]; names(point1) = c("x","y")
				point2 = segments[i,c("x2","y2")]; names(point2) = c("x","y")
				points = rbind(point1,point2); linesList = list()
				linesList[[1]] = Lines(list(Line(points)),1)
				spatialLine = SpatialLines(linesList)
				crs(spatialLine) = crs(shapefiles[[1]])
				for (j in 1:length(shapefiles))
					{
						crossingBarrierEvent = 0
						intersections = gIntersection(shapefiles[[j]], spatialLine)
						crossingBarrierEvent = length(intersections)
						if (!is.null(intersections)) crossingBarrierEvent = 1
						crossingBarrierEvents[i,j] = crossingBarrierEvent
						if ((showingPlots == TRUE)&(crossingBarrierEvent>0))
							{
								lines(spatialLine, col="red", lwd=0.5)
							}
						if ((showingPlots == TRUE)&(crossingBarrierEvent==0))
							{
								lines(spatialLine, col="green3", lwd=0.2)
							}
					}
			}
		if (showingPlots == TRUE) dev.off()
		return(crossingBarrierEvents)
	}
Ns_obs_list = list(); Ns_ran_list = list()
for (h in 1:length(datasets))
	{
		data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
		data2 = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), head=T)
		segments_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_segments.csv"))
		
		# 3.3.1. Determining the number of edges crossing a barrier
		
		shapefiles = list(international_borders)
		crossingBarriers_obs = countingCrossingBarrierEvents(segments_obs, shapefiles)
		write.csv(crossingBarriers_obs, paste0("Outbreak_results/",datasets[h],"_dir/Obs_crossings.csv"), row.names=F, quote=F)
		registerDoMC(10); buffer = list()
		buffer = foreach(i = 1:nberOfRandomisations) %dopar% {
		# for (i in 1:nberOfRandomisations) {
				segments_ran = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_segments_",i,".csv"), header=T)
				crossingBarriers_ran = countingCrossingBarrierEvents(segments_ran, shapefiles)
				fileName = paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_crossings_",i,".csv")
				write.csv(crossingBarriers_ran, fileName, row.names=F, quote=F)
				i
			}
		
		# 3.3.2. Statistical test to assess the impact of barriers
		
		Ns_ran = list()
		crossingBarriers_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_crossings.csv"), head=T)
		counters = rep(0, dim(crossingBarriers_obs)[2])
		observed_Ns = matrix(nrow=1, ncol=dim(crossingBarriers_obs)[2])
		randomised_Ns = matrix(nrow=nberOfRandomisations, ncol=dim(crossingBarriers_obs)[2])
		colnames(observed_Ns) = colnames(crossingBarriers_obs)
		colnames(randomised_Ns) = colnames(crossingBarriers_obs)
		for (i in 1:nberOfRandomisations)
			{
				crossingBarriers_ran = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_crossings_",i,".csv"), header=T)
				for (j in 1:dim(crossingBarriers_obs)[2])
					{
						N_obs = sum(crossingBarriers_obs[,j]); observed_Ns[1,j] = N_obs
						N_ran = sum(crossingBarriers_ran[,j]); randomised_Ns[i,j] = N_ran
						if (N_ran <= N_obs) counters[j] = counters[j] + 1
					}
			}
		pValue = t(counters/nberOfRandomisations)
		cat(datasets[[h]]," - International_borders, p-value = ",pValue,"\n",sep="")
		saveRDS(observed_Ns, paste0("Outbreak_results/",datasets[h],"_dir/N_obs_values.rds"))
		saveRDS(randomised_Ns, paste0("Outbreak_results/",datasets[h],"_dir/N_ran_values.rds"))
		Ns_obs_list[[h]] = observed_Ns; Ns_ran_list[[h]] = randomised_Ns
		if (showingPlots == TRUE)
			{
				par(oma=c(0,0,0,0), mar=c(5,5,4,2))
				plot(density(randomised_Ns), ylab="density", xlab="N", main=datasets[h], cex.main=1, cex.lab=1)
				abline(v=observed_Ns, col="red")
			}
	}

	# 3.4. Analysis of the impact of local immunity on the dispersal location
	
sumOfImmunityAtEndLocations = function(segments, shapefile, immunityValues, startingDate)
	{
		Is = rep(NA, dim(segments)[1]); pointsWithoutAdmin2 = c(); savingTemporaryFile = FALSE
		for (i in 1:dim(segments)[1])
			{
				pt.x = segments[i,"x2"]; pt.y = segments[i,"y2"]; indices1 = c()
				for (j in 1:length(shapefile@polygons))
					{
						for (k in 1:length(shapefile@polygons[[j]]@Polygons))
							{
								pol.x = shapefile@polygons[[j]]@Polygons[[k]]@coords[,1]
								pol.y = shapefile@polygons[[j]]@Polygons[[k]]@coords[,2]
								if (point.in.polygon(pt.x,pt.y,pol.x,pol.y) == 1)
									{
										indices1 = c(indices1, j)
									}
							}
					}
				if (length(indices1) != 1)
					{
						cat(paste0("    ",i,", length of indices1 = ",length(indices1)),"\n",sep="")
						if (length(indices1) == 0)
							{
								pointsWithoutAdmin2 = rbind(pointsWithoutAdmin2, cbind(pt.x, pt.y))
							}
					}	else	{
						admin2_code = shapefile@data[indices1,"ADM2_CODE"]
						date = decimal_date(date_decimal(startingDate)+days(segments[i,"dt"]))
						indices2 = which((immunityValues[,"admin2_code"]==admin2_code)&(date>=immunityValues[,"startMonth"])&(date<immunityValues[,"endMonth"]))
						if (length(indices2) != 1)
							{
								cat(paste0("    ",i,", length of indices2 = ",length(indices2)),"\n",sep="")
							}	else	{
								Is[i] = immunityValues[indices2,"immunity"]
							}
					}
			}
		if ((savingTemporaryFile)&(!is.null(pointsWithoutAdmin2)))
			{
				pdf(paste0("TEMP.pdf"), width=5, height=2.3)
				par(mfrow=c(1,1), mar=c(0.0,0.0,0.0,0.4), oma=c(0.4,0.0,0.6,1.6), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				plot(polyExtent, ann=F, axes=F, col=NA, border=NA)
				plot(popDensity_008_log, col=colorRampPalette(brewer.pal(9,"Greys"))(17)[2:10], box=F, axes=F, legend=F, add=T)
				plot(international_polygons, add=T, lwd=0.2, col=NA, border="gray30", lty=1)
				points(pointsWithoutAdmin2, cex=0.3, pch=16, col="red"); points(pointsWithoutAdmin2, cex=0.3, pch=1, lwd=0.3, col="gray30")
				rect(xmin(template), ymin(template), xmax(template), ymax(template), xpd=T, lwd=0.2, border="gray30")
				dev.off()
			}
		if (!is.null(pointsWithoutAdmin2)) cat("    Warning: ",dim(pointsWithoutAdmin2)[1]," ending points without an admin-2 match\n",sep="")
		I = mean(Is[!is.na(Is[])])
		return(I)
	}
Is_obs_list = list(); Is_ran_list = list()
for (h in 1:length(datasets))
	{
		data1 = read.csv(paste0("Outbreak_results/",datasets[h],".csv"), head=T)
		data2 = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Data_KDE_95.csv"), head=T)
		segments_obs = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/Obs_segments.csv"))
		
		# 3.4.1. Computing the sum of immunity values at the end location of all segments
		
		shapefile = shapefile("Africa_adm2_shp/Africa_admin2.shp")
		immunityValues = read.csv("SimImmunity.csv", head=T) # @Darlan: to be edited with your CSV file containing real immunity estimates
		startMonths = as.POSIXct(paste0(immunityValues[,"date"],"-01 00:00")); endMonths = startMonths+months(1)
		startEndMonths = cbind(decimal_date(startMonths), decimal_date(endMonths))
		colnames(startEndMonths) = c("startMonth","endMonth"); immunityValues = cbind(immunityValues, startEndMonths)
		startingDate = data1[which(data1[,"collection_date"]==min(data1[,"collection_date"])),"collection_date"]
		I_obs = sumOfImmunityAtEndLocations(segments_obs, shapefile, immunityValues, startingDate)
		registerDoMC(10); Is_ran = rep(NA, nberOfRandomisations); buffer = list()
		buffer = foreach(i = 1:nberOfRandomisations) %dopar% {
		# for (i in 1:nberOfRandomisations) {
				segments_ran = read.csv(paste0("Outbreak_results/",datasets[h],"_dir/All_ran_segments/Ran_segments_",i,".csv"), header=T)
				sumOfImmunityAtEndLocations_ran = sumOfImmunityAtEndLocations(segments_ran, shapefile, immunityValues, startingDate)
				sumOfImmunityAtEndLocations_ran
			}
		for (i in 1:length(buffer)) Is_ran[i] = buffer[[i]]
		
		# 3.4.2. Statistical test to assess the impact of immunity values at ending locations
		
		counter = 0
		for (i in 1:nberOfRandomisations)
			{
				if (I_obs <= Is_ran[i]) counter = counter + 1
			}
		pValue = counter/nberOfRandomisations
		cat(datasets[[h]]," - Impact of immunity (I stat test), p-value = ",pValue,"\n",sep="")
		Is_obs_list[[h]] = I_obs; Is_ran_list[[h]] = Is_ran
		if (showingPlots == TRUE)
			{
				par(oma=c(0,0,0,0), mar=c(5,5,4,2))
				plot(density(Is_ran), ylab="density", xlab="N", main=datasets[h], cex.main=1, cex.lab=1)
				abline(v=I_obs, col="red")
			}
	}

