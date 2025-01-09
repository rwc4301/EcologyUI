library(jsonlite)
library(httpgd)
library(unigd)
library(ggplot2)
library(plumber)

#* @apiTitle R Plot API

#* Enable CORS
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "*")
  res$setHeader("Access-Control-Allow-Headers", "*")
  res$setHeader("Access-Control-Allow-Credentials", "true")
  # res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  # res$setHeader("Access-Control-Allow-Headers", "Content-Type")

  if (req$REQUEST_METHOD == "OPTIONS") {
    res$status <- 200
    return(list())
  }

  plumber::forward()
}

#* Get a hello message
#* @get /hello
function() {
  list(message = "R Server Connected")
}

#* Generate a histogram plot
#* @param bins The number of bins for the histogram
#* @post /plot
function(bins = 30) {
  # Initialize httpgd if not already running
  if (!dev.cur() == 2) hgd()
  
  # Create the plot
  p <- ggplot(faithful, aes(x = waiting)) +
    geom_histogram(bins = as.numeric(bins), fill = "#007bc2", color = "white") +
    labs(x = "Waiting time to next eruption (in mins)",
         y = "Count",
         title = "Histogram of waiting times")
  
  # Render to httpgd device
  print(p)
  
  # Return base64 encoded plot
  list(plotData = ugd_render())
}

