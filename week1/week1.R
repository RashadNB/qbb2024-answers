library(ggplot2)

# Load coverage data from text file
coverage <- read.delim("~/qbb2024-answers/week1/coverage3.txt", header = FALSE)
poissonCoverage = coverage %>% dplyr::mutate(Poisson = 1000000*dpois(0:999999,3))

# Create a sequence of integer x-values for Poisson (discrete distribution)
lambda <- 30
x_values <- 0:max(coverage$V1)  # Use range of coverage data
poisson_probs <- dpois(x_values, lambda)

maxcoverage = max(coverage$V1)

normal_estimates = 1000000*dnorm(0:maxcoverage, 30, sqrt(30))

# Create a data frame for Poisson distribution
poisson_df <- data.frame(x = x_values, y = poisson_probs)

# Plot the histogram of the coverage data
ggplot() + 
  # Histogram from the coverage data
  geom_histogram(data = coverage, 
                 mapping = aes(x = V1, y = after_stat(count)), 
                 binwidth = 1, 
                 fill = "skyblue", 
                 color = "white") +
  
  # Add labels and theme
  labs(title = "Coverage Distribution with Poisson(λ = 3) and Normal Distribution", 
       x = "Coverage", 
       y = "Counts") +
  theme_minimal() +
  
  # Overlay Poisson distribution
  geom_smooth(data = poisson_df, aes(x = x, y = y * length(coverage$V1) / sum(poisson_probs), color = "Poisson")) +
  geom_smooth(aes(x = 0:maxcoverage, y = normal_estimates, color = "Normal")) +
  scale_color_manual(values = c(Poisson = "blue" , Normal = "red")) +
  ylab("Frequency")
  

              