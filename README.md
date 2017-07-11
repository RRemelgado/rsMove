## rsMove
Bridging Remote Sensing and Movement Ecology with R.

<br>

### 1. Why Develop rsMove?

<p align="justify">
In the scope of movement ecology, Global Positioning Systems (GPS) have evolved significantely offering an unique insight into the animal behavior. But understanding this behavior is dependent on our ability to compreeend the underlying environmental conditions that guides it. In this context, remote sensing becomes a fundamental tool. It provides information on the spatial and temporal variability of the landscape and provides us the means to understand the impact of environmental change over animal behavior. However, linking remote sensing and animal movement can be tricky due to the differences in the spatial and temporal scales at which they are acquired (Figure 1). As a consequence, a simple point-to-raster query becomes insuficient creating a demand for data handling methods that are sensitive to the constraints imposed by remote sensing. <i>rsMove</i> Answers to this issue providing tools to query and analyze movement data using remote sensing that are sensitive to the spatial and temporal limitations of satellite based environmental predictors.
</p>

<br>

<p align="center">
  <img width="566" height="291" src="http://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs40462-015-0036-7/MediaObjects/40462_2015_36_Fig1_HTML.gif">

  Figure 1 - Scale differences between animal movement and remotely-sensed data (<a href="https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-015-0036-7">Neuman et al, 2015</a>)  
</p>

<br>

### 2. Some generic tools (but mostly not)
<p align="justify">
The development of packages such as <i>raster</i> and <i>sp</i> opened the door for the use of remote sensing within R. they provide generic tools to process spatial data as well as an efficient approaches to handle the large datasets that are characteristic of the field of remote sensing. As a result, <i>rsMove</i> aims not to replicate the work done within these packages but rather extend its applicability to the particular issues that characterize its usage within movement ecology. In this section, we provide some examples on the main applicabilities of this package.
</p>

#### 2.1. Example I - Finding hotspots to find test sites

#### 2.2. Example II - Sampling in time
<p align="justify">
Wen using multi-temporal remote sensing data to understand animal movement patterns, querying remote sensing data can be challenging due to the dynamic spatial and temporal nature of animal tracking data. The function <i>dataQuery()</i> offers an interface to link remote sensing and movement data by selecting the most adequate temporal information for each GPS record. The function offers the choice for an exact match of an approximate match. While the first option only allows the selection of remote sensing data collected at the same time as the GPS record, the second is more permissive searching for the nearest clear pixel within a user defined temporal buffer.
</p>

```R
moveData <- shapefile(system.file('extdata', 'konstanz_20130805-20130811.shp', package="rsMove"))
rsStack <- stack(list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=T))
rsQuery <- dataQuery(xy=moveData,img=rsStack)
```
<p align="justify">
Independently of our ability to select temporal information, remote sensing suffers from an old and inescapable evil: cloud cover.

NOTE: TALK ABOUT INTERPOLATION

</p>

#### 2.3. Example III - Sampling on the move

### 3. Instalation
This gitHub is used as a basis for the improvement of *rsMove*. A stable release is avalible on CRAN and can installed with

```R
install.packages('rsMove')
```


