################################################
# Install the libraries
install.packages("trend")
install.packages("TSA")
install.packages("forecast")

##################################################
# Load the libraries
library(trend)
library(TSA)
library(forecast)

##################################################
# Import CSV file
file_path <- file.choose()  # Select the CSV file
dat <- read.csv(file_path, header = TRUE, row.names = 1)
da <- ts(dat)

##################################################
# Handle missing values
da <- na.omit(da)

##################################################
# MAN-KENDALL Test
apply_mk_test <- function(data) {
  for (i in 1:ncol(data)) {
    cat("MAN-KENDALL Test for column:", colnames(data)[i], "\n")
    print(mk.test(data[, i]))
  }
}
apply_mk_test(da)

##################################################
# COX-STUART Test
apply_cs_test <- function(data) {
  for (i in 1:ncol(data)) {
    cat("COX-STUART Test for column:", colnames(data)[i], "\n")
    
    # Check if there is variation in the data
    if (length(unique(data[, i])) == 1) {
      cat("Error: Entire sample contains identical values\n")
    } else {
      print(cs.test(data[, i]))
    }
  }
}
# Apply COX-STUART Test
apply_cs_test(da)

#################################################
# Plot Partial Autocorrelation AND Autocorrelation
for (i in 1:ncol(da)) {
  # some columns may have 0, we can skip them
  if (sd(da[, i]) != 0) {
    # Partial Autocorrelation
    pacf(da[, i], main = paste("Partial Autocorrelation of", colnames(da)[i]))
    
    # Autocorrelation
    acf(da[, i], main = paste("Autocorrelation of", colnames(da)[i]))
  } else {
    cat("Skipping plotting for column:", colnames(da)[i], "as it contains all zeroes.\n")
  }
}

#################################################
# Calculate the Periodicity
# Initialize an empty list to store results
results_list <- list()

# Loop over each column
for (i in 1:ncol(da)) {
  # Calculate periodogram
  periodogram_result <- periodogram(da[, i], main = colnames(da)[i])
  
  # Extract relevant data
  dd <- data.frame(freq = periodogram_result$freq, spec = periodogram_result$spec)
  order <- dd[order(-dd$spec), ]
  top <- head(order, 5)
  time <- 1 / top$f
  
  # Store results in the list
  results_list[[i]] <- list(top, time)
}

for (j in 1:ncol(da)) {
  cat("Results for column", colnames(da)[j], ":\n")
  print(results_list[[j]][[1]])  # Top frequencies
  print(results_list[[j]][[2]])  # Corresponding periods
}

#################################################
# Seasonal decomposition using LOESS method
# Initialize an empty list to store seasonal decomposition results
decomposition_list <- list()

# Loop over each column
for (i in 1:ncol(da)) {
  # Create time series object
  ts_data <- ts(da[, i], frequency = 12, start = c(2004, 1))
  
  # Seasonal decomposition
  stl_result <- stl(ts_data, s.window = "periodic")
  
  # Store decomposition results in the list
  decomposition_list[[i]] <- stl_result
}

# Plot seasonal decomposition for each column
for (j in 1:ncol(da)) {
  # Extract decomposition results
  stl_result <- decomposition_list[[j]]
  
  # Plot
  plot(stl_result, main = paste("Seasonal decomposition in", colnames(da)[j]))
}

###############################################
# Seasonal decomposition using TBATS method
# Initialize an empty list to store TBATS results
tbats_results <- list()

# Loop over each column
for (i in 1:ncol(da)) {
  # Create time series object
  ts_data <- ts(da[, i], frequency = 12, start = c(2010, 1))
  
  # Fit TBATS model
  tbats_model <- tbats(ts_data, use.box.cox = NULL, use.trend = NULL, 
                       use.damped.trend = NULL, seasonal.periods = c(3, 4, 6), 
                       use.arma.errors = TRUE, num.cores = 2, bc.lower = 0, bc.upper = 1)
  
  # Store TBATS results
  tbats_results[[i]] <- tbats_model
  
  # Plot seasonal decomposition
  plot(tbats_model, main = paste("Seasonal decomposition in", colnames(da)[i]))
}

##############################################

