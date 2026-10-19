# This program adjusts inputs for a mean and SD in the model to work in a truncated normal distribution.
# Author: Mingshu Huang (jenniferhuangmsh@gmail.com or for urgent questions, Mingshu.Huang@carelon.com)
# Co-Author: Virginia Talbot (vtalbot@mgh.harvard.edu)
# App can be run at this website: https://virginia-apps.shinyapps.io/AgeFitter/
# Original version: https://mingshuhuang.shinyapps.io/test/

library(shiny)
library(msm)

# this is the function that builds the web interface
ui <- fluidPage(
  titlePanel("Age Fitter"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "mean_input", label = "Mean", value = 480, min = 0, max = 1200), 
      numericInput(inputId = "mean_min", label = "Lower Bound (minimum age)", value = 0, min = 0),
      numericInput(inputId = "mean_max", label = "Upper Bound (maximum age)", value = 1200, min = 1200),
      numericInput(inputId = "sd_input", label = "SD", value = 120, min = 0),
      numericInput(inputId = "weight", label = "Enter weight for mean", value = 1, min = 0),
      textOutput(outputId = "warning"),
      actionButton(inputId = "calculate", label = "Calculate!")
    ),
    mainPanel(
      plotOutput(outputId = "hist"),
      verbatimTextOutput(outputId = "summary"),
      tableOutput(outputId = "stat")
    )
  )
)

# this is a function that calculated the mean and sd of the truncated normal given the data
optimization <- function(mu, sd, lower, upper, weight) {
  
  # define the left and right bounder of the grid for mu
  mu_left <- max(lower, mu - 5)
  mu_right <- min(upper, mu + 5)
  # define the left and right bounder of the grid for sd
  sd_left <- sd
  sd_right <- sd * 4
  # this is the weight given to the mean when minimizing the difference between the desired and fitted mean and sd values
  mu_weight <- weight
  
  # create a sequence of small incrementation within the boundary
  mu_seq <- seq(from = mu_left, to = mu_right, by = 0.1)
  sd_seq <- seq(from = sd_left, to = sd_right, by = 0.1)
  # create the mu and sd matrix for later computation
  mu_len <- length(mu_seq)
  sd_len <- length(sd_seq)
  mu_matrix <- matrix(rep(mu_seq, sd_len), nrow = mu_len, ncol = sd_len)
  sd_matrix <- t(matrix(rep(sd_seq, mu_len), nrow = sd_len, ncol = mu_len))
  # calculated the truncated mean and sd for each value of mean and sd in the grid
  left_z <- (lower - mu_matrix) / sd_matrix
  right_z <- (upper - mu_matrix) / sd_matrix
  left_d <- dnorm(left_z)
  right_d <- dnorm(right_z)
  left_p <- pnorm(left_z)
  right_p <- pnorm(right_z)
  trunc_mu <- mu_matrix - sd_matrix * (left_d - right_d) / (left_p - right_p)
  trunc_sd <- sd_matrix * sqrt(1 - (left_z * left_d - right_z * right_d) / (left_p - right_p) - ((left_d - right_d) / (left_p - right_p))^2)
  # calculate the difference between the desired values and the computed values
  weighted_diff <- mu_weight * abs(mu - trunc_mu) + abs(sd - trunc_sd)
  # return the normal parameters and the fitted truncated mean and sd
  position <- which(weighted_diff==min(weighted_diff), arr.ind=T)
  # zoom in for more accurate estimates
  new_mu <-  mu_matrix[position]
  new_sd <-  sd_matrix[position]
  new_mu_left <- max(mu_left, new_mu - 0.1)
  new_mu_right <- min(mu_right, new_mu + 0.1)
  new_sd_left <- max(sd_left, new_sd - 0.1)
  new_sd_right <- min(sd_right, new_sd + 0.1)
  new_mu_seq <- seq(from = new_mu_left, to = new_mu_right, by = 0.01)
  new_sd_seq <- seq(from = new_sd_left, to = new_sd_right, by = 0.01)
  new_mu_len <- length(new_mu_seq)
  new_sd_len <- length(new_sd_seq)
  new_mu_matrix <- matrix(rep(new_mu_seq, new_mu_len), nrow = new_mu_len, ncol = new_sd_len)
  new_sd_matrix <- t(matrix(rep(new_sd_seq, new_mu_len), nrow = new_sd_len, ncol = new_mu_len))
  new_left_z <- (lower - new_mu_matrix) / new_sd_matrix
  new_right_z <- (upper - new_mu_matrix) / new_sd_matrix
  new_left_d <- dnorm(new_left_z)
  new_right_d <- dnorm(new_right_z)
  new_left_p <- pnorm(new_left_z)
  new_right_p <- pnorm(new_right_z)
  new_trunc_mu <- new_mu_matrix - new_sd_matrix * (new_left_d - new_right_d) / (new_left_p - new_right_p)
  new_trunc_sd <- new_sd_matrix * sqrt(1 - (new_left_z * new_left_d - new_right_z * new_right_d) / (new_left_p - new_right_p) - ((new_left_d - new_right_d) / (new_left_p - new_right_p))^2)
  new_weighted_diff <- mu_weight * abs(mu - new_trunc_mu) + abs(sd - new_trunc_sd)
  new_position <- which(new_weighted_diff==min(new_weighted_diff), arr.ind=T) 
  return(c(new_norm_mu=new_mu_matrix[new_position], new_norm_sd=new_sd_matrix[new_position], new_tn_mu=new_trunc_mu[new_position], new_tn_sd=new_trunc_sd[new_position], norm_mu=mu_matrix[position], norm_sd=sd_matrix[position], tn_mu=trunc_mu[position], tn_sd=trunc_sd[position]))
}


# this is the server side function that takes in user's web input, runs it, and gives results back to the user
server <- function(input, output) { 
  observeEvent(input$calculate, {
    if (input$mean_min < 0) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The age lower bound must be a nonnegative number."
      })
      output$warning <- renderText({})
    }
    else if (input$mean_input <= 0) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The mean age must be a positive number."
      })
      output$warning <- renderText({})
    }
    else if (input$mean_max > 1200) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The age upper bound must be no greater than 1200."
      })
      output$warning <- renderText({})
    }
    else if (input$mean_min >= input$mean_max) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The age upper bound must be greater that the age lower bound."
      })
      output$warning <- renderText({})
    }
    else if (input$mean_input <= input$mean_min) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The mean age must be greater than the age lower bound."
      })
      output$warning <- renderText({})
    }    
    else if (input$mean_input >= input$mean_max) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The mean age must be smaller than the age upper bound."
      })
      output$warning <- renderText({})
    }
    else if (input$sd_input <= 0) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The SD must be a positive number."
      })
      output$warning <- renderText({})
    }
    else if (input$weight <= 0) {
      output$hist <- renderPlot({})
      output$stat <- renderTable({})
      output$summary <- renderText({
        "Error: The weight must be a positive number."
      })
      output$warning <- renderText({})
    }
    else {
      result <- optimization(input$mean_input, input$sd_input, input$mean_min, input$mean_max, input$weight)
      target_mean <- input$mean_input
      target_sd <- input$sd_input
      lb <-input$mean_min
      ub <- input$mean_max
      weight <- input$weight
      # below is to check the fitted values by simulating 100000 random values from the resulted truncated normal distribution and calculate the mean and sd
      data <- rtnorm(100000, mean = result[1], sd = result[2], lower=lb, upper=ub)
      summary_value <- round(c(target_mean, mean(data), target_sd, sd(data), result[1], result[2]), digit=2)
      summary_name <- c("targeted mean", "fitted mean", "targeted SD", "fitted SD", "Mean to enter in CEPAC", "SD to enter in CEPAC")
      smry <- data.frame(summary_name, summary_value)
      colnames(smry) <- c("name", "value")
      output$hist <- renderPlot({
        title <- paste("Result of 100,000 simulation of truncated normal\nlower bound = ", lb, ", upper bound = ", ub, ", weight = ", weight, sep="")
        hist(data, main = title, xlim = c(lb, ub), xlab = "Age", breaks = 20)
      })
      output$stat <- renderTable({
        smry
      })
      output$summary <- renderPrint({
        summary(data, digits=4)
      })
      if (input$weight > 10)
      {
        output$warning <- renderText({
          "Warning: you are giving more than 10 times the weight to mean as to SD."
        })
      }
      else if (input$weight < 0.1)
      {
        output$warning <- renderText({
          "Warning: you are giving less than 1/10 the weight to mean as to SD."
        })
      }
      else
      {
        output$warning <- renderText({})
      }
    }
  })
}

# this says: run what you have built above!
shinyApp(ui = ui, server = server)
