## rsMove
Bridging Remote Sensing and Movement Ecology with R.

<br>

### 1. Why Develop rsMove?

<p align="justify">
In the scope of movement ecology, Global Positioning Systems (GPS) have evolved significantely offering an unique insight into the animal behavior. But understanding this behavior is dependent on our ability to compreeend the underlying environmental conditions that guides it. In this context, remote sensing becomes a fundamental tool. It provides information on the spatial and temporal variability of the landscape and provides us the means to understand the impact of environmental change over animal behavior. However, linking remote sensing and animal movement can be troublesome due to the differences in the spatial and temporal scales at which they are acquired (Figure 1). As a consequence, methods that are sensitive to the constraints imposed by remote sensing in the analysis of animal movement are required. <i>rsMove</i> Answers to this issue providing tools to query and analyze movement data using remote sensing.
</p>

<br>

<p align="center">
  <img width="566" height="291" src="http://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs40462-015-0036-7/MediaObjects/40462_2015_36_Fig1_HTML.gif">
</p>

<p align="center">
   Figure 1 - Scale differences between animal movement and remotely-sensed data (<a href="https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-015-0036-7">Neuman et al, 2015</a>)  
</p>

<br>

### 2. Instalation
This gitHub is used as a basis for the improvement of *rsMove*. A stable release is avalible on CRAN and can installed with

```R
install.packages('rsMove')
```

### 3. Example - Finding hotspots to find test sites
<p align="justify">
Within this section, we provide an example on the combine use of the functions  <i>sampleMove()</i> and  <i>hotMove()</i> for the identification of areas of interest that can serve as test sites. For this purpose, we rely on White Stork movement data (DOI: <a href="10.5441/001/1.78152p3q">10.5441/001/1.78152p3q</a>) which was colected by the Max Planck Institute for Ornithologie (MPIo) and is accessible through <a href="https://www.movebank.org/">MoveBank</a>. We focused on data collected between June and December of 2013 which refers to the first migration of 13 juveniles between Radofzell, in Germany, and the Gibraltar narrow, at the coast of Spain (Figure 2).
</p>

<br>

 <p align="center"><img width="600" height="290" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_2.jpg"></p>
 
<p align="center">Figure 2 - Movement track from 13 Juvenile white Storks between Germany and the Gribraltar narrow.</p>

```R

# read movement data
shp <- shapefile("/WhiteStork_Konstanz.shp")

# Extract time information and convert to date format
ot = as.Date(shp@data$date)

# Derive unique identifiers
ui <- unique(shp@data$ind)

# sample for each unique identifier
os = vector('list', length(ui))
for (i in 1:length(ui)) {
    ind = which(shp@data$ind==ui[i])
    os[[i]] <- sampleMove(xy=shp[ind,], ot=as.Date(shp@data$dare[ind]), error=7, method="deg")
}

```

<p align="center"><img width="600" height="290" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_3.jpg"></p>
 
<p align="center">Figure 2 - Comparison between original samples (in black) and the samples selected by <i>sampleMove</i> (in red).</p>

```R

# derive shapefile from samples of each unique ID
x = unlist(sapply(os, function(x){x@coords[,1]}))
y = unlist(sapply(os, function(x){x@coords[,2]}))
shp <- SpatialPoints(cbind(x,y), proj4string=crs(shp))

# identify sample regions
hm <- hotMove(xy=s, pxr=0.2, shp=T)

```

<p align="center"><img width="600" height="290" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_4.jpg"></p>
 
<p align="center">Figure 2 - Polygons of different colors represent unique sample regions identified with <i>hotMove</i> (in red).</p>

