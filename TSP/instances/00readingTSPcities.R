# http://www.math.uwaterloo.ca/tsp/world/countries.html
# http://www.math.uwaterloo.ca/tsp/world/

#setwd on data folder!


download.TSP <- function(sample){
  url.base <- "http://www.math.uwaterloo.ca/tsp/world/"
  data <- paste0(url.base, sample, ".tsp")
  tour <- paste0(url.base, sample, ".tour")
  file.data <- paste0(sample, ".tsp")
  file.tour <- paste0(sample, ".tour")
  download.file(data, file.data)
  download.file(tour, file.tour)
}

download.TSP("lu980")
download.TSP("mo14185")
download.TSP("nu3496")
download.TSP("mu1979")
download.TSP("pm8079")
download.TSP("qa194")
download.TSP("rw1621")
download.TSP("sw24978")
download.TSP("tz6117")
download.TSP("uy734")
download.TSP("vm22775")
download.TSP("wi29")
download.TSP("ym7663")
download.TSP("zi929")
download.TSP("ca4663")
download.TSP("it16862")
download.TSP("ar9152")
download.TSP("bm33708")
download.TSP("ch71009")
download.TSP("ei8246")
download.TSP("it16862")
download.TSP("kz9976")

#http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/

