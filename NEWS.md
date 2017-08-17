=================================================================================
## rsMove 0.2.2
=================================================================================

### New:
  * Introduced new functions
    * sMoveRes()
    * tMoveRes()
    * specVar()
    * getEnv()
    * proSat()
    * moveCloud()
    * satTime()
    * moveReduce()
    * plotMove()
 * imgInt(), moveSeg(), timeDir(), spaceDir(), linInt() now accept data frames as inputs.
 * hotMoveStats() now reports on each temporal segment found within each region and returns the corresponding sample indices.
 * spaceDir() accepts categorical data.
 * spaceDir() and moveSeg() use an optional spatial buffer.
 * moveSeg(), timeDir(), spaceDir() and hotMoveStats() return plots.
 * dataQuery() requires an adjustable buffer specifying the size of search window before and after the target date.

### Changes
  * spaceDirSample() was renamed as spaceDir() to avoid confusion regarding its purpose.
  * timeDirSample() was renamed as timeDir() to avoid confusion regarding its purpose.
  * moveMovel() now requires data frames for presences and absences instead of shapefiles.
  * moveSeg(), spaceDir() and timeDir() no longer samples to single pixels internaly.
 
### Fixes:
 * imgInt() allows that the interpolation is limited to the records in a point shapefile.
 * moseSeg() now recognizes categorical data.
 * BackSample() now assings background samples to the correct pixels.

=================================================================================
## rsMove 0.1.0
=================================================================================

Initial release to CRAN (2017-07-15) with the following functions:
 * backSample()
 * dataQuery()
 * hotMove()
 * hotMoveStats()
 * imgInt()
 * labelSample()
 * modelApply()
 * moveModel()
 * moveSeg()
 * poly2sample()
 * rsComposite()
 * sampleMove()
 * segRaster()
 * spaceDirSample()
 * timeDirSample()
