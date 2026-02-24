# This file suppresses NOTEs about "no visible binding for global variable"
# during R CMD check.

utils::globalVariables(c(
  # Column names / Tidy evaluation variables
  "Attribute.Name", "imp", "mod", "dfC", "ev", "staticVars", "bathyR",

  # Operators from other packages
  ".", "%>%", "%do%", "%dopar%"
))
