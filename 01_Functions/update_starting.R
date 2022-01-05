update_starting <- function(parameters, starting, ag, bg) {

  ag <- parameters$ag_biomass_min +
    (parameters$ag_biomass_max - parameters$ag_biomass_min) * ag
  
  bg <- parameters$bg_biomass_min +
    (parameters$bg_biomass_max - parameters$bg_biomass_min) * bg
  
  starting$ag_biomass <- ag
  
  starting$bg_biomass <- bg
  
  return(starting)
}
