library(shiny)
library(jsonlite)
library(httpgd)
library(ggplot2)

# Add this before the server function
addResourcePath("api", ".")

# Custom HTTP response handler for CORS
corsHandler <- function(handler) {
  function(req) {
    if (req$REQUEST_METHOD == "OPTIONS") {
      # Handle preflight requests
      list(
        status = 200L,
        headers = list(
          "Access-Control-Allow-Origin" = "*",
          "Access-Control-Allow-Methods" = "GET, POST, OPTIONS",
          "Access-Control-Allow-Headers" = "Content-Type",
          "Access-Control-Allow-Credentials" = "true"
        ),
        body = ""
      )
    } else {
      # Handle actual requests
      res <- handler(req)
      res$headers[["Access-Control-Allow-Origin"]] <- "*"
      res$headers[["Access-Control-Allow-Methods"]] <- "GET, POST, OPTIONS"
      res$headers[["Access-Control-Allow-Headers"]] <- "Content-Type"
      res$headers[["Access-Control-Allow-Credentials"]] <- "true"
      res
    }
  }
}

server <- function(input, output, session) {
  # Initialize httpgd
  hgd()

  # session$onSessionEnded(function() {
  #   hgd_destroy()
  # })

  # Hello World endpoint using observeEvent
  observeEvent(input$getHello, {
    response <- list(message = "Hello from R!")
    session$sendCustomMessage("helloResponse", response)
    output$text <- renderText({
        "Button was clicked!"
    })
    print("Hello from R!")
  })
  
  # Plot endpoint using observeEvent
  observeEvent(input$getPlot, {
    req(input$bins)  # Ensure bins parameter exists
    
    # Create example ggplot
    p <- ggplot(faithful, aes(x = waiting)) +
      geom_histogram(bins = input$bins, fill = "#007bc2", color = "white") +
      labs(x = "Waiting time to next eruption (in mins)",
           y = "Count",
           title = "Histogram of waiting times")
    
    # Render plot to httpgd device
    print(p)
    
    # Get plot as base64 PNG and send response
    plot_data <- hgd_base64()
    session$sendCustomMessage("plotResponse", list(plotData = plot_data))
  })
  
  # Wrap session$registerDataObj with CORS handler
  # session$registerDataObj(
  #   name = "plot",
  #   data = NULL,
  #   filterFunc = corsHandler(function(data, req) {
  #     # Parse the request body
  #     body <- rawToChar(req$rook.input$read())
  #     params <- fromJSON(body)
      
  #     # Create example ggplot
  #     p <- ggplot(faithful, aes(x = waiting)) +
  #       geom_histogram(bins = params$bins, fill = "#007bc2", color = "white") +
  #       labs(x = "Waiting time to next eruption (in mins)",
  #            y = "Count",
  #            title = "Histogram of waiting times")
      
  #     # Render plot to httpgd device
  #     print(p)
      
  #     # Get plot as base64 PNG
  #     plot_data <- hgd_base64()
      
  #     # Return JSON response
  #     list(
  #       status = 200L,
  #       headers = list("Content-Type" = "application/json"),
  #       body = toJSON(list(plotData = plot_data))
  #     )
  #   })
}

ui <- fluidPage(
  tags$div(
    h1("Shiny API Example"),
  ),
  actionButton("getHello", "Get Hello Action"),
  tags$button(id = "getHello", "Get Hello"),
  tags$button(id = "getPlot", "Get Plot"),
  textOutput("text"),
  tags$script(HTML('
    // Handle CORS preflight
    Shiny.addCustomMessageHandler("helloResponse", function(message) {
      window.dispatchEvent(new CustomEvent("helloResponse", { detail: message }));
    });
    
    Shiny.addCustomMessageHandler("plotResponse", function(message) {
      window.dispatchEvent(new CustomEvent("plotResponse", { detail: message }));
    });
  '))
)

shinyApp(ui = ui, server = server)