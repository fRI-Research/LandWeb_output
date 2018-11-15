defineModule(sim, list(
  name = "LandWeb_output",
  description = "Summarize the output for the LandWeb natural range of variation (NRV).",
  keywords = c("LandWeb", "NRV"),
  authors = person("Yong", "Luo", email = "yong.luo@canada.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("1.3.2"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandWeb_output.Rmd"),
  reqdPkgs = list("data.table", "raster", "SpaDES.tools"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("summaryInterval", "numeric", 50, NA, NA, "This describes summary interval for this module"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "summaryPeriod", objectClass = "numeric",
                 desc = "a numeric vector contains the start year and end year of summary",
                 sourceURL = NA),
    expectsInput(objectName = "vegLeadingProportion", objectClass = "numeric",
                 desc = "a number that define whether a species is lead for a given pixel",
                 sourceURL = NA),
    expectsInput(objectName = "cohortData", objectClass = "data.table",
                 desc = "age cohort-biomass table hooked to pixel group map by pixelGroupIndex at
                 succession time step, this is imported from forest succession module",
                 sourceURL = NA),
    expectsInput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step, this is imported from
                 forest succession module",
                 sourceURL = NA),
    expectsInput("species", "data.table", "Columns: species, speciesCode, Indicating several features about species",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("vegTypeMapGenerator", "function", "converts a species table and cohortdata and pixelGroupMap to vegTypeMap raster")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "vegTypeMap", objectClass = "Raster", desc = NA)
  )
))

doEvent.LandWeb_output <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {

    sim <- scheduleEvent(sim, 0, "LandWeb_output", "initialConditions", eventPriority = 1)
    sim <- scheduleEvent(sim, 0, "LandWeb_output", "allEvents", eventPriority = 7.5)
    sim <- scheduleEvent(sim, sim$summaryPeriod[1], "LandWeb_output", "allEvents",
                         eventPriority = 7.5)
  } else if (eventType == "initialConditions") {
    plotVTM(sim$specieslayers, vegLeadingProportion = sim$vegLeadingProportion,
            speciesEquivalency = sim$speciesEquivalency)
  } else if (eventType == "allEvents") {
    if (time(sim) >= sim$summaryPeriod[1] &
        time(sim) <= sim$summaryPeriod[2]) {
      sim <- AllEvents(sim)
      sim <- scheduleEvent(sim,  time(sim) + P(sim)$summaryInterval,
                           "LandWeb_output", "allEvents", eventPriority = 7.5)
    }
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.


plotVTM <- function(speciesStack = NULL, vtm = NULL, vegLeadingProportion, speciesEquivalency) {
  if (is.null(vtm)) {
    if (!is.null(speciesStack))
      vtm <- Cache(pemisc::makeVegTypeMap, speciesStack, vegLeadingProportion)
    else
      stop("plotVTM requires either a speciesStack of percent cover or a vegetation type map (vtm)")
  }

  vtmTypes <- factorValues(vtm, seq(minValue(vtm), maxValue(vtm)), att = "Species")[[1]]
  setColors(vtm, vtmTypes) <-
    equivalentName(vtmTypes, df = speciesEquivalency, "cols")
  #facVals <- factorValues(vtm, vtm[], att = "Species")[[1]]
  facVals <- pemisc::factorValues2(vtm, vtm[], att = "Species", na.rm = TRUE)
  df <- data.table(species = as.character(facVals), stringsAsFactors = FALSE)
  df <- df[!is.na(df$species)]
  df$species <- equivalentName(df$species, speciesEquivalency, "shortNames")
  df$cols <- equivalentName(df$species, speciesEquivalency, "cols")

  cols2 <- df$cols
  names(cols2) <- df$species
  initialLeadingPlot <- ggplot(data = df, aes(species, fill = species)) +
    scale_fill_manual(values=cols2) +
    geom_bar(position = 'stack') +
    #labs(x = "Year", y = "Biomass by species") +
    theme(legend.text = element_text(size = 6), legend.title = element_blank(),
          axis.text = element_text(size = 6))


  Plot(initialLeadingPlot, title = c("Initial leading types"))
  Plot(vtm, title = "Initial leading types")

}
### template for your event1
AllEvents <- function(sim) {
  sim$vegTypeMap <- sim$vegTypeMapGenerator(sim$species, sim$cohortData, sim$pixelGroupMap,
                                            sim$vegLeadingProportion)

  # vegetation type summary
  # if(is.null(sim$LandMine$vegTypeMapGenerator)) { # This may be produced in a specific fire module
  #   species <- sim$species
  #   species[species == "Pinu_sp" | species == "Pinu_sp", speciesGroup := "PINU"]
  #   species[species == "Betu_pap" | species == "Popu_bal"|
  #             species == "Popu_tre" | species == "Lari_lar", speciesGroup := "DECI"]
  #   species[species == "Pice_mar" | species == "Pice_gla", speciesGroup := "PICE"]
  #   cohortdata <- sim$cohortData
  #   shortcohortdata <- setkey(cohortdata, speciesCode)[setkey(species[,.(speciesCode, speciesGroup)],
  #                                                             speciesCode), nomatch = 0]
  #   shortcohortdata[, totalB := sum(B, na.rm = TRUE), by = pixelGroup]
  #   shortcohortdata <- shortcohortdata[,.(speciesGroupB = sum(B, na.rm = TRUE),
  #                                         totalB = mean(totalB, na.rm = TRUE)),
  #                                      by = c("pixelGroup", "speciesGroup")]
  #   shortcohortdata[,speciesProportion:=speciesGroupB/totalB]
  #   shortcohortdata[speciesGroup == "PINU" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 1]# pine leading
  #   shortcohortdata[speciesGroup == "DECI" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 2]# deciduous leading
  #   shortcohortdata[speciesGroup == "PICE" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 3]# spruce leading
  #   shortcohortdata[is.na(speciesLeading), speciesLeading := 0]
  #   shortcohortdata[,speciesLeading:=max(speciesLeading, na.rm = TRUE), by = pixelGroup]
  #   shortcohortdata <- unique(shortcohortdata[,.(pixelGroup, speciesLeading)], by = "pixelGroup")
  #   shortcohortdata[speciesLeading == 0, speciesLeading := 4] # 4 is mixed forests
  #   attritable <- data.table(ID = unique(shortcohortdata$speciesLeading))
  #   attritable[ID == 1, Factor := "Pine leading"]
  #   attritable[ID == 2, Factor := "Deciduous leading"]
  #   attritable[ID == 3, Factor := "Spruce leading"]
  #   attritable[ID == 4, Factor := "Mixed"]
  #   pixelGroupMap <- sim$pixelGroupMap
  #   vegTypeMap <- rasterizeReduced(shortcohortdata, pixelGroupMap, "speciesLeading")
  #   vegTypeMap <- setValues(vegTypeMap, as.integer(getValues(vegTypeMap)))
  #   levels(vegTypeMap) <- as.data.frame(attritable)
  #   projection(vegTypeMap) <- projection(sim$pixelGroupMap)
  #   sim$vegTypeMap <- vegTypeMap
  # } else {
  #   sim$vegTypeMap <- sim$LandMine$vegTypeMapGenerator(sim$species, sim$cohortData, sim$pixelGroupMap,
  #                                             sim$vegLeadingProportion)
  # }
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  if (!suppliedElsewhere("summaryPeriod", sim))
    sim$summaryPeriod <- c(1000, 1500)
  if (!suppliedElsewhere("vegTypeMapGenerator", sim)) { # otherwise created in LandMine
    sim$vegTypeMapGenerator <- function(species, cohortdata, pixelGroupMap, vegLeadingProportion) {
      species[species == "Pinu_ban" | species == "Pinu_con" | species == "Pinu_sp", speciesGroup := "PINU"]
      species[species == "Betu_pap" | species == "Popu_bal" | species == "Popu_tre" |
                species == "Lari_lar", speciesGroup := "DECI"]
      species[species == "Pice_mar" , speciesGroup := "PICE_MAR"]
      species[species == "Pice_gla", speciesGroup := "PICE_GLA"]
      species[species == "Abie_sp" , speciesGroup := "ABIE"]
      #cohortdata <- sim$cohortData
      shortcohortdata <- setkey(cohortdata, speciesCode)[setkey(species[, .(speciesCode, speciesGroup)],
                                                                speciesCode), nomatch = 0]
      shortcohortdata[, totalB := sum(B, na.rm = TRUE), by = pixelGroup]
      shortcohortdata <- shortcohortdata[, .(speciesGroupB = sum(B, na.rm = TRUE),
                                             totalB = mean(totalB, na.rm = TRUE)),
                                         by = c("pixelGroup", "speciesGroup")]
      shortcohortdata[,speciesProportion := speciesGroupB/totalB]

      speciesLeading <- NULL
      Factor <- NULL
      ID <- NULL
      pixelGroup <- NULL
      speciesProportion <- NULL
      speciesGroup <- NULL
      speciesCode <- NULL
      totalB <- NULL
      B <- NULL
      speciesGroupB <- NULL

      shortcohortdata[speciesGroup == "PINU" & speciesProportion > vegLeadingProportion,
                      speciesLeading := 1]# pine leading
      shortcohortdata[speciesGroup == "DECI" & speciesProportion > vegLeadingProportion,
                      speciesLeading := 2]# deciduous leading
      shortcohortdata[speciesGroup == "PICE_MAR" & speciesProportion > vegLeadingProportion,
                      speciesLeading := 3]# spruce leading
      shortcohortdata[speciesGroup == "PICE_GLA" & speciesProportion > vegLeadingProportion,
                      speciesLeading := 4]# spruce leading
      shortcohortdata[is.na(speciesLeading), speciesLeading := 0]
      shortcohortdata[, speciesLeading := max(speciesLeading, na.rm = TRUE), by = pixelGroup]
      shortcohortdata <- unique(shortcohortdata[, .(pixelGroup, speciesLeading)], by = "pixelGroup")
      shortcohortdata[speciesLeading == 0, speciesLeading := 5] # 5 is mixed forests
      attritable <- data.table(ID = sort(unique(shortcohortdata$speciesLeading)))
      attritable[ID == 1, Factor := "Pine leading"]
      attritable[ID == 2, Factor := "Deciduous leading"]
      attritable[ID == 3, Factor := "Black spruce leading"]
      attritable[ID == 4, Factor := "White spruce leading"]
      attritable[ID == 5, Factor := "Mixed"]
      vegTypeMap <- rasterizeReduced(shortcohortdata, pixelGroupMap, "speciesLeading", "pixelGroup")
      vegTypeMap <- setValues(vegTypeMap, as.integer(getValues(vegTypeMap)))
      levels(vegTypeMap) <- as.data.frame(attritable)
      projection(vegTypeMap) <- projection(pixelGroupMap)
      vegTypeMap
    }
  }
  if (!suppliedElsewhere("vegLeadingProportion", sim)) {
    sim$vegLeadingProportion <- 0.8
  }
  if (!suppliedElsewhere("cohortData", sim))
    sim$cohortData <- data.table()
  if (!suppliedElsewhere("pixelGroupMap", sim))
    sim$pixelGroupMap <- raster()
  if (!suppliedElsewhere("species", sim)) {
    localSpeciesFilename <- file.path(dataPath(sim), "speciesTraits.csv")
    if (!file.exists(localSpeciesFilename)) {
      mm <- moduleMetadata(currentModule(sim), getPaths()$modulePath)$inputObjects
      download.file(subset(mm, objectName == "species")$sourceURL,
                    destfile = localSpeciesFilename)
    }
    sim$species <- read.csv(localSpeciesFilename, header = TRUE,
                            stringsAsFactors = FALSE) %>%
      data.table()
  }
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
