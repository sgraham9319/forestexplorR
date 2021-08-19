#' Abiotic data for 15 forest stands.
#' 
#' A data frame containing measurements of five abiotic variables in each of
#' 15 forest stands located in Mount Rainier National Park, WA, USA. These
#' forest stands are part of the Pacific Northwest Permanent Sample Plot Network 
#' (\url{http://pnwpsp.forestry.oregonstate.edu/})
#' 
#' @format A data frame with 15 rows and 6 variables:
#' \describe{
#'   \item{stand_id}{unique identification code of each stand}
#'   \item{precip_mm}{annual precipitation in mm}
#'   \item{temp_C}{mean annual temperature in degrees Celsius}
#'   \item{elev_m}{stand elevation in m}
#'   \item{aet_mm}{actual evapotranspiration in mm}
#'   \item{pet_mm}{potential evapotranspiration in mm}
#' }
#' @source Precipitation and temperature data were obtained from PRISM 
#' (\url{https://prism.oregonstate.edu/}) and represent the 1981-2010 30yr
#' normals. Evapotranspiration data were derived from precipitation and
#' temperature data using the method outlined in Ford et al. 2017. Can. J. For.
#' Res. 47, 53-62. Elevation data were obtained from the Pacific Northwest
#' Permanent Sample Plot Network
#' (\url{http://pnwpsp.forestry.oregonstate.edu/data})
"stand_abiotic"
#' 