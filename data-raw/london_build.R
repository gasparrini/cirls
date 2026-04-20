# Read raw data
data <- read.csv("data-raw/london.csv")

# Add some variables for ease of use in examples
london <- dplyr::mutate(data,
  date = as.Date(date, format = "%d/%m/%Y"),
  doy = as.numeric(format(date, "%j")),
  dow = weekdays(date))

# Rearrange columns
london <- dplyr::relocate(london, doy, dow, .after = date)

# Save
save(london, file = "data/london.rda")
