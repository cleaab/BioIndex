#' MEDITS coordinates in decimal degrees
#'
#' @param Data ddata frame of TA table
#'
#' @return the function return the same dataframe with the coordinates converted in the decimal degrees format
#' @export convert_coordinates

convert_coordinates<-function(Data)  {

  lat_start=Data$SHOOTING_LATITUDE
  lon_start= Data$SHOOTING_LONGITUDE
  lat_end=Data$HAULING_LATITUDE
  lon_end= Data$HAULING_LONGITUDE
  LatStartDeg = floor(floor(lat_start)/100);
  LonStartDeg = floor(floor(lon_start)/100);
  LatStartMin=(lat_start-LatStartDeg*100)/60
  LonStartMin=(lon_start-LonStartDeg*100)/60
  LatEndDeg = floor(floor(lat_end)/100);
  LonEndDeg = floor(floor(lon_end)/100);
  LatEndMin=(lat_end-LatEndDeg*100)/60
  LonEndMin=(lon_end-LonEndDeg*100)/60

  lat_start2= LatStartDeg + LatStartMin
  lon_start2 = LonStartDeg + LonStartMin
  lat_end2 = LatEndDeg + LatEndMin
  lon_end2 = LonEndDeg + LonEndMin
  Data$lat_start = lat_start2
  Data$lon_start = lon_start2
  Data$lat_end = lat_end2
  Data$lon_end = lon_end2
  return(Data)
}
