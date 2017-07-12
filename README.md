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

### 3. Example - Selecting test sites through hotspot detection
<p align="justify">
Within this section, we provide an example on the combine use of the functions  <i>sampleMove()</i>,  <i>hotMove()</i> and <i>hotMoveStats()</i> for the identification of areas of interest that can serve as test sites. For this purpose, we rely on White Stork movement data (DOI: <a href="10.5441/001/1.78152p3q">10.5441/001/1.78152p3q</a>) which was colected by the Max Planck Institute for Ornithologie (MPIo) and is accessible through <a href="https://www.movebank.org/">MoveBank</a>. We focused on data collected between June and December of 2013 which refers to the first migration of 13 juveniles between Radofzell, in Germany, and the Gibraltar narrow, at the coast of Spain (Figure 2). Due to size restrictions, this example data is not provided in the CRAN version of the package. However, you can obtain it <a href="">within this gitHub</a>. This is pre-processed version of the original data provided through MoveBank.
</p>

<br>

<p align="center"><img width="605" height="315" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_2.jpg"></p>
 
<p align="center">Figure 2 - Movement track from 13 Juvenile white Storks between Germany and the Gribraltar narrow.</p>

<br>

This data was collected with a temporal resolution of 5 minutes and an acquisition error of 3.6 meters. Due to its fine resolution, it provides valuable information on the behavior of this species allowing us to understand how it interacts with its environment during different behavior stages. Figure 3 shows the movement of the White Storks analyzed within this example while they were nesting in Germany. This GIF, derived with the <a href="https://github.com/cran/moveVis">moveVis</a> package, shows that the animals spent most of their time in small density urban areas where their nests were located. But we also see how they move towards the outskirts of the urban fabric and into agricultural land when in search for food. As a result, this data gives us an insight into the differences in resource selection by the White Stork when it's resting and when it's feeding. But it also includes a lot of random information related to the travels between the nesting and feeding sites. Additionally, we can see some exploratory flights where the animals don't seem to focus on a particular area. While this info can be relevant to understand how the animal moves, from a remote sensing perspective, these samples account for noise that hinders our ability to accurately map relevant resources.

<br>

<p align="center"><img width="605" height="315" src="https://github.com/RRemelgado/rsMove/blob/master/ws_animation.jpg"></p>
 
<p align="center">Figure 3 - Animation of White Stork movement tracks (Radofzell, Germany) narrow.</p>

<br>

<p align="center"><img width="605" height="315" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_2.jpg"></p>
 
<p align="center">Figure 2 - Movement track from 13 Juvenile white Storks between Germany and the Gribraltar narrow.</p>

<br>

As a consequence, before we can use this data to model the suitability of environmental resources, we first need to indetify samples that represent areas of interest and exclude random movements. <i>sampleMove()</i> offers a simple appproach to achieve this goal. Looking at continuous movements, the function evaluates the distance between consecutive points to identify the start of a stop. If the distance between consecutive points is lesser than the predefined error (in this case 7.2 meters, which is two times the error of the GPS sensor) the first points is set as a pointer. Then, the following observations are compared against this pointer and added to the segment if the distance between them is lower than the error. If not, the segment is closed and the function provides an average for the coordinates and timestamps as well as an estimate of the elapsed time within it. Below you can see how to read in the data and use <i>sampleMove()</i as well as the spatial distribution of the retrieved samples (Figure 3).


```R

# read movement data
shp <- shapefile("./WhiteStork_Germany.shp")

# Extract time information and convert to date format
ot = as.Date(shp@data$date)

# Derive unique identifiers
ui <- unique(shp@data$ind)

# sample for each unique identifier
os = vector('list', length(ui))
for (i in 1:length(ui)) {
    ind = which(shp@data$ind==ui[i])
    os[[i]] <- sampleMove(xy=shp[ind,], ot=as.Date(shp@data$dare[ind]), error=7.2, method="deg")
}

```

<p align="center"><img width="605" height="315" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_3.jpg"></p>
 
<p align="center">Figure 3 - Comparison between original samples (in black) and the samples selected by <i>sampleMove</i> (in red).</p>

```R

# derive shapefile from samples of each unique ID
x = unlist(sapply(os, function(x){x@coords[,1]}))
y = unlist(sapply(os, function(x){x@coords[,2]}))
shp <- SpatialPoints(cbind(x,y), proj4string=crs(shp))

# identify sample regions
hm <- hotMove(xy=s, pxr=0.2, shp=T)

```

<p align="center"><img width="605" height="315" src="https://github.com/RRemelgado/rsMove/blob/master/Figure_4.jpg"></p>
 
<p align="center">Figure 4 - Polygons of different colors represent unique sample regions identified with <i>hotMove</i> (in red).</p>

