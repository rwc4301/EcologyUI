library(jsonlite)
library(httpgd)
library(unigd)
library(ggplot2)
library(plumber)

# Initialize httpgd device
hgd()
hgd_url()

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

#* Get device state
#* @get /state
function() {
  state <- ugd_state()
  details <- hgd_details()
  url <- hgd_url()

  list(
    num_plots = state$hsize,
    update_id = state$upid,
    active = state$active,
    client = state$client,
    host = details$host,
    port = details$port,
    token = details$token,
    status = details$status
  )
}

#* Clear the device
#* @post /clear
function() {
  dev.off()
  hgd()
  list(success = TRUE)
}

#* Get a hello message
#* @get /hello
function() {
  list(message = "R Server Connected")
}

#* Get a hello message
#* @get /do
function() {
  print(ggplot(data = mtcars, aes(x = mpg, y = wt)) + geom_point())
  list(message = "R Server Connected")

}

#* Generate a histogram plot
#* @param bins The number of bins for the histogram
#* @post /plot
function(bins = 30) {
  # # Initialize httpgd if not already running
  # if (!dev.cur() == 2) hgd()
  
  # Create the plot
  p <- ggplot(faithful, aes(x = waiting)) +
    geom_histogram(bins = as.numeric(bins), fill = "#007bc2", color = "white") +
    labs(x = "Waiting time to next eruption (in mins)",
         y = "Count",
         title = "Histogram of waiting times")
  
  # Render to httpgd device
  print(p)

  # Return the current state
  list(success = TRUE)
  
  # # Return base64 encoded plot
  # list(plotData = ugd_render())
}

