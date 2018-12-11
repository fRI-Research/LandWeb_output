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
  reqdPkgs = list("data.table", "raster", "SpaDES.tools", "PredictiveEcology/pemisc"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("sppEquivCol", "character", "LandWeb", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("summaryInterval", "numeric", 50, NA, NA, "This describes summary interval for this module"),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    desc = "a number that define whether a species is leading for a given pixel"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput("cohortData", "data.table",
                 desc = paste("age cohort-biomass table hooked to pixel group map by pixelGroupIndex at",
                              "succession time step, this is imported from forest succession module"),
                 sourceURL = ""),
    expectsInput("species", "data.table",
                 desc = "Columns: species, speciesCode, Indicating several features about species",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "updated community map at each succession time step",
                 sourceURL = ""),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See pemisc::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon to use as the study area,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput("summaryPeriod", "numeric",
                 desc = "a numeric vector contains the start year and end year of summary",
                 sourceURL = "")
  ),
  outputObjects = bind_rows(
    createsOutput("vegTypeMap", "Raster", desc = NA)
  )
))

doEvent.LandWeb_output <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandWeb_output", "initialConditions",
                         eventPriority = 1)
    sim <- scheduleEvent(sim, 0, "LandWeb_output", "allEvents", eventPriority = 7.5)
    sim <- scheduleEvent(sim, sim$summaryPeriod[1], "LandWeb_output", "allEvents",
                         eventPriority = 7.5)
  } else if (eventType == "initialConditions") {
    plotVTM(speciesStack = stack(raster::mask(sim$speciesLayers, sim$rasterToMatch)),
            vegLeadingProportion = P(sim)$vegLeadingProportion,
            sppEquiv = sim$sppEquiv,
            sppEquivCol = P(sim)$sppEquivCol,
            title = "Initial Types")
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
#   - keep event functions short and clean, modularize by calling subroutines from section below.

AllEvents <- function(sim) {
  sim$vegTypeMap <- vegTypeMapGenerator(sim$species, sim$cohortData, sim$pixelGroupMap,
                                        P(sim)$vegLeadingProportion)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects", "function:spades")
  cPath <- cachePath(sim)
  dPath <- asPath(dataPath(sim), 1)

  if (!suppliedElsewhere("summaryPeriod", sim))
    sim$summaryPeriod <- c(1000, 1500)

  if (!suppliedElsewhere("cohortData", sim))
    sim$cohortData <- data.table()

  if (!suppliedElsewhere("pixelGroupMap", sim))
    sim$pixelGroupMap <- raster()

  if (!suppliedElsewhere("species", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath, cacheTags)
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    data("sppEquivalencies_CA", package = "pemisc", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)

    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## add default colors for species used in model
    defaultCols <- RColorBrewer::brewer.pal(6, "Accent")
    LandRNames <- c("Pice_mar", "Pice_gla", "Popu_tre", "Pinu_sp", "Abie_sp")
    sim$sppEquiv[LandR %in% LandRNames, cols := defaultCols[-4]]
    sim$sppEquiv[EN_generic_full == "Mixed", cols := defaultCols[4]]
  }

  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    speciesLayersList <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyAreaLarge,
                               speciesList = sim$speciesList,
                               # thresh = 10,
                               url = extractURL("speciesLayers"),
                               cachePath = cachePath(sim),
                               userTags = c(cacheTags, "speciesLayers"))

    #options(opts)
    writeRaster(speciesLayersList$speciesLayers,
                file.path(outputPath(sim), "speciesLayers.grd"),
                overwrite = TRUE)
    sim$speciesLayers <- speciesLayersList$speciesLayers
    #sim$speciesList <- speciesLayersList$speciesList ## not used in this module
  }
  return(invisible(sim))
}
