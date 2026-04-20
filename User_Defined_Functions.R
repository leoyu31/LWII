# Header ------------------------------------------------------------------
# Author: Liyao Yu
# Date: 2025-Feb-20
# User defined functions


# Theme ---------------------------------------------------------------------------------------
theme_Yu <- theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.2),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6)
  )

# Return the abbreviated species name -------------------------------------
abbr_species <- function(species_name) {
  vapply(species_name, function(name) {
    words <- unlist(strsplit(name, " "))
    if (length(words) < 2) {
      stop(paste0(words, " is invalid species name. It should have at least two words."))
    }
    paste0(substr(words[1], 1, 1), ". ", words[2])
  }, FUN.VALUE = character(1))  # Ensures output is character
}

# Convert DoY to Month ----------------------------------------------------
doy_to_month <- function(doy, year) {
  date <- as.Date(doy - 1, origin = paste0(year, "-01-01"))
  month <- format(date, "%m")
  return(as.numeric(month))
}

# Calculate how many days there are in a month ----------------------------
days_in_month <- function(year, month) {
  if (any(month < 1 | month > 12)) {
    stop("Invalid input: 'month' must be an integer between 1 and 12.")
  }
  # Adjust for December (month = 12)
  next_month_year <- ifelse(month == 12, year + 1, year)
  next_month <- ifelse(month == 12, 1, month + 1)
  # Calculate the last day of the current month
  last_day <- as.Date(paste(next_month_year, next_month, "1", sep = "-"), 
                      format = "%Y-%m-%d") - 1
  return(as.numeric(format(last_day, "%d")))
}

# Exclude outliers with a given quantile (> 50%) ----------------------------
exclude_outliers_vec <- function(x, y) {
  x <- as.numeric(x)  # ensure numeric
  if (y <= 0.5)  y <- 1 - y
  lower_threshold <- quantile(x, 1 - y, na.rm = TRUE)
  upper_threshold <- quantile(x, y, na.rm = TRUE)
  x[x < lower_threshold | x > upper_threshold] <- NA
  return(x)
}
# Function for data frame
exclude_outliers_df_str <- function(df, colname, y) {
  df[[colname]] <- exclude_outliers_vec(df[[colname]], y)
  return(df)
}

# Calculate mean, sd and cv at a rolling base -----------------
roll_stats <- function(x, width = 5) {
  m <- rollapply(x, width, mean, na.rm = TRUE, align = "right", fill = NA)
  s <- rollapply(x, width, sd,   na.rm = TRUE, align = "right", fill = NA)
  tibble(mean = m, sd = s, cv = s / m)
}