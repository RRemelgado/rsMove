## rsMove
Bridging Remote Sensing and Movement Ecology with R.

<br>

### Why develop rsMove?

<p align="justify">
In the scope of movement ecology, Global Positioning Systems (GPS) have evolved significantly offering a unique insight into the animal behavior. But understanding this behavior is dependent on our ability to comprehend the underlying environmental conditions that guides it. In this context, remote sensing becomes a fundamental tool. It provides information on the spatial and temporal variability of the landscape and provides us the means to understand the impact of environmental change over animal behavior. However, linking remote sensing and animal movement can be troublesome due to the differences in the spatial and temporal scales at which they are acquired (Figure 1). As a consequence, methods that are sensitive to the constraints imposed by remote sensing in the analysis of animal movement are required. <i>rsMove</i> Answers to this issue providing tools to query and analyze movement data using remote sensing.
</p>

<br>

<p align="center">
  <img width="566" height="291" src="http://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs40462-015-0036-7/MediaObjects/40462_2015_36_Fig1_HTML.gif">
</p>

<p align="center">
<sub>Figure 1 - Scale differences between animal movement and remotely-sensed data (<a href="https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-015-0036-7">Neuman et al, 2015</a>)</sub>
</p>

<br>

### Installation
This gitHub is used as a basis for the improvement of *rsMove*. A stable release is available on CRAN and can installed with:

```R
install.packages('rsMove')
```
### Examples
* <a href="https://github.com/RRemelgado/README_data/blob/master/rsMove/example_1.md">Finding the hotspots! Ojective study site selection.</a>


<br>

### Aknowledgements
<p align="justify">
This initiative is part of the <a href="http://www.fernerkundung.geographie.uni-wuerzburg.de/en/lehrstuehle_und_arbeitsgruppen/department_of_remote_sensing/research/projects/current_projects/opt4environment//">Opt4Environment</a> project and was funded by the German Aerospace Center (DLR) on behalf of the Federal Ministry for Economic Affairs and Energy (BMWi) with the research grant <b<50 EE 1403</b>.For news on other developments at the Department of Remote Sensing of the University of Würzburg, clike <a href="http://remote-sensing.eu/">here</a>
</p>

### Other packages developed by our department
* <a href="http://bleutner.github.io/RStoolbox/">RStookbox</a>
* <a href="https://github.com/cran/moveVis/">moveVis</a>
