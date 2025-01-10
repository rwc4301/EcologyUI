library(fs)
library(plumber)

run_api <- function() {
  pr <- plumb("plumber.R")
  pr$run(port = 3839)
}

watch_files <- function(path, interval = 2) {
  print("Watching for changes in files...")
  last_time <- dir_info(path)$modification_time
  print(last_time)
  repeat {
    Sys.sleep(interval)
    new_time <- dir_info(path)$modification_time
    if (any(new_time != last_time)) {
      print("File changed!")
      run_api()
      last_time <- new_time
    }
  }
}

run_api()
# watch_files(".")