## rsMove
Bridging Remote Sensing and Movement Ecology with R.

<br>

### 1. Why Develop rsMove?
In the scope of movement ecology, Global Positioning Systems (GPS) have evolved significantely offering an unique insight into the animal behavior. But understanding this behavior is dependent on our ability to compreeend the underlying environmental conditions that guides it. In this context, remote sensing becomes a fundamental tool. It provides information on the spatial and temporal variability of the landscape and provides us the means to understand the impact of environmental change over animal behavior. However, linking remote sensing and animal movement can be tricky due to the differences in the spatial and temporal scales at which they are acquired (Figure 1). As a consequence, a simple point-to-raster query becomes insuficient creating a demand for data handling methods that are sensitive to the constraints imposed by remote sensing. rsMove Answers to this issue providing tools to query and analyze movement data using remote sensing that are sensitive to the spatial and temporal limitations of satellite based environmental predictors.

<br>

<p align="center">
  <img width="566" height="291" src="http://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs40462-015-0036-7/MediaObjects/40462_2015_36_Fig1_HTML.gif">
</p>

Figure 1 - Scale difference between animal movement and remotely-sensed data ([Neuman et al, 2015](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-015-0036-7))

<br>

### 2. Some generic tools (but mostly not)
The development of packages such as *raster* and *sp* open the door for the use of remote sensing within R. they provide generic tools to process spatial data as well as an efficient approaches to handle the large datasets that are characteristic of the field of remote sensing. As a result, *rsMove* aims not to replicate the work done within these packages but rather extend its applicability to the particular issues that characterize its usage within movement ecology. In this section, we discuss some of the main applicabilities of this package.

#### 2.1. Spatial and temporal querying
Wen using multi-temporal remote sensing data to understand animal movement patterns, querying remote sensing data can be challenging due to the dynamic nature of animal tracking data. The function *dataQuery()* offers an interface to link remote sensing and movement data by selecting the most adequate temporal information for each GPS record. The function offers the choice for an exact match of an approximate match. While the first option only allows the selection of remote sensing data collected at the same time as the GPS record, the second is more permissive searching for the nearest clear pixel within a user defined temporal buffer.

<br>

```R
moveData <- shapefile(system.file('extdata', 'konstanz_20130805-20130811.shp', package="rsMove"))
rsStack <- stack(list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=T))
rsQuery <- dataQuery(xy=moveData,img=rsStack)
```




Thus explaining its behavior is dependent on our ability to represent the observed changes in the environment.


While animal tracking data is acquired on very fine temporal scales (e.g. minutes, hours) satellite data is limited by its overpass time with open-access satellites such as the Moderate Resolution Spectrodadiometer (MODIS) offering daily acquisitions. Additionaly, more than often, the quality quality of satellite data is hindered by clouds that lead to annoying temporal gaps.

The function *dataQuery()* extends on the *extract()* function provided by the *raster* package and 
