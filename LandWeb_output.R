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
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput("cohortData", "data.table",
                 desc = paste("age cohort-biomass table hooked to pixel group map by pixelGroupIndex at",
                              "succession time step, this is imported from forest succession module"),
                 sourceURL = ""),
    expectsInput("fireReturnInterval","Raster",
                 desc = "A raster layer that is a factor raster, with at least 1 column called fireReturnInterval, representing the fire return interval in years"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "updated community map at each succession time step",
                 sourceURL = ""),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "this raster contains two pieces of information: Full study area with fire return interval attribute", ## TODO: is this correct?
                 sourceURL = NA),
    expectsInput("rstTimeSinceFire", "Raster", "a time since fire raster layer", NA),
    expectsInput("species", "data.table",
                 desc = "Columns: species, speciesCode, Indicating several features about species",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("sppColors", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in sim$speciesEquivalency[[sim$sppEquivCol]],",
                              "and should also contain a color for 'Mixed'"),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See pemisc::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("standAgeMap", "RasterLayer",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureStandVolume.tar"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon to use as the study area,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA),
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
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandWeb_output", "otherPlots",
                         eventPriority = 1)
    sim <- scheduleEvent(sim, 0, "LandWeb_output", "allEvents", eventPriority = 7.5)
    sim <- scheduleEvent(sim, sim$summaryPeriod[1], "LandWeb_output", "allEvents",
                         eventPriority = 7.5)
  } else if (eventType == "initialConditions") {
    plotVTM(speciesStack = raster::mask(sim$speciesLayers, sim$studyAreaReporting) %>% stack(),
            vegLeadingProportion = P(sim)$vegLeadingProportion,
            sppEquiv = sim$sppEquiv,
            sppEquivCol = P(sim)$sppEquivCol,
            colors = sim$sppColors,
            title = "Initial Types")

    ## plot initial age map
    ageMap <- raster::mask(sim$standAgeMap, sim$studyAreaReporting) %>% stack()
    Plot(ageMap, title = "Initial stand ages")
  } else if (eventType == "allEvents") {
    if (time(sim) >= sim$summaryPeriod[1] &
        time(sim) <= sim$summaryPeriod[2]) {
      sim <- AllEvents(sim)
      sim <- scheduleEvent(sim,  time(sim) + P(sim)$summaryInterval,
                           "LandWeb_output", "allEvents", eventPriority = 7.5)
    }
  } else if (eventType == "otherPlots") {
    ## average age by FRI polygon
    tsfMap <- raster::mask(sim$rstTimeSinceFire, sim$studyAreaReporting)
    fris <- unique(na.omit(sim$fireReturnInterval[]))
    names(fris) <- fris
    tsfs <- vapply(unname(fris), function(x) {
      ids <- which(sim$fireReturnInterval[] == x)
      unname(mean(tsfMap[ids], na.rm = TRUE))
    }, numeric(1))
    polys <- sim$fireReturnInterval
    tsfDF <- data.frame(time = as.numeric(times(sim)$current),
                        meanAge = unname(tsfs),
                        FRI = as.factor(unname(fris)))
    mod$tsfOverTime <- rbind(mod$tsfOverTime, tsfDF)
    mod$tsfOverTime <- mod$tsfOverTime[!is.na(mod$tsfOverTime$meanAge),]

    if (length(unique(mod$tsfOverTime$time)) > 1) {
      gg_tsfOverTime <- ggplot(mod$tsfOverTime,
                               aes(x = time, y = meanAge, col = FRI, ymin = 0)) +
        geom_line(size = 1.5) +
        theme(legend.text = element_text(size = 14))

      firstPlot <- identical(time(sim), P(sim)$.plotInitialTime + P(sim)$.plotInterval)
      title1 <- if (firstPlot) "Average age (TSF) by FRI polygon" else ""

      Plot(gg_tsfOverTime, title = title1, addTo = "ageOverTime")
    }

    ## schedule future plots
    sim <- scheduleEvent(sim, times(sim)$current + P(sim)$.plotInterval, "LandWeb_output",
                         "otherPlots", eventPriority = 1)
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

AllEvents <- function(sim) {
  sim$vegTypeMap <- vegTypeMapGenerator(sim$cohortData, sim$pixelGroupMap,
                                        P(sim)$vegLeadingProportion,
                                        colors = sim$sppColors)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyArea <- randomStudyArea(seed = 1234)
  }

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaLarge <- sim$studyArea
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaLarge <- sim$studyArea
  }

  if (!suppliedElsewhere("fireReturnInterval", sim))
    stop("fireReturnInterval map must be supplied.")

  if (!suppliedElsewhere("rasterToMatch", sim))
    stop("rasterToMatch must be supplied.")

  if (!suppliedElsewhere("rasterToMatchReporting")) {
    sim$rasterToMatchReporting <- sim$rasterToMatch
  }

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
    if (!is.null(sim$sppColors))
      stop("If you provide sppColors, you MUST also provide sppEquiv")
    sim$sppColors <- pemisc::sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                       newVals = "Mixed", palette = "Accent")
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

  if (!suppliedElsewhere("standAgeMap", sim)) {
    sim$standAgeMap <- Cache(prepInputs, #notOlderThan = Sys.time(),
                             targetFile = basename(standAgeMapFilename),
                             archive = asPath(c("kNN-StructureStandVolume.tar",
                                                "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.zip")),
                             destinationPath = dPath,
                             url = extractURL("standAgeMap"),
                             fun = "raster::raster",
                             studyArea = sim$studyAreaLarge,
                             rasterToMatch = sim$rasterToMatch,
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = TRUE, overwrite = TRUE,
                             userTags = c("stable", currentModule(sim)))
  }

  return(invisible(sim))
}
