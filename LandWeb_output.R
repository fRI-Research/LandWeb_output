defineModule(sim, list(
  name = "LandWeb_output",
  description = "Summarize the output for the LandWeb natural range of variation (NRV).",
  keywords = c("LandWeb", "NRV"),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Yong", "Luo", email = "yluo1@lakeheadu.ca", role = "aut"),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(LandWeb_output = numeric_version("1.3.3")),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandWeb_output.Rmd"),
  reqdPkgs = list("data.table", "raster", "SpaDES.tools",
                  "PredictiveEcology/LandR@development (>= 1.0.9.9001)",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    defineParameter("mixedType", "numeric", 2,
                    desc = paste("How to define mixed stands: 1 for any species admixture;",
                                 "2 for deciduous > conifer. See ?vegTypeMapGenerator.")),
    defineParameter("sppEquivCol", "character", "LandWeb", NA, NA,
                    desc = "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("summaryInterval", "numeric", 50, NA, NA,
                    desc = "This describes summary interval for this module"),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    desc = "a number that define whether a species is leading for a given pixel"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".plots", "character", default = "object",
                    desc = paste("Passed to `types` in `Plots` (see `?Plots`). There are a few plots that are made within this module, if set.",
                                 "Note that plots (or their data) saving will ONLY occur at `end(sim)`.",
                                 "If `NA`, plotting is turned off completely (this includes plot saving).")),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput("cohortData", "data.table",
                 desc = paste("age cohort-biomass table hooked to pixel group map by `pixelGroupIndex` at",
                              "succession time step, this is imported from forest succession module."),
                 sourceURL = NA),
    expectsInput("fireReturnInterval", "Raster",
                 desc = paste("A raster layer that is a factor raster,",
                              "with at least 1 column called fireReturnInterval,",
                              "representing the fire return interval in years.")),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "updated community map at each succession time step",
                 sourceURL = NA),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = paste("this raster contains two pieces of information:",
                              "Full study area with fire return interval attribute."), ## TODO: is this correct?
                 sourceURL = NA),
    expectsInput("rstTimeSinceFire", "Raster",
                 desc = "a time since fire raster layer",
                 sourceURL = NA),
    expectsInput("speciesTable", "data.table",
                 desc = "Columns: species, speciesCode, Indicating several features about species",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in `sim$speciesEquivalency[[sim$sppEquivCol]]`,",
                              "and should also contain a color for 'Mixed'"),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See `LandR::sppEquivalencies_CA`.",
                 sourceURL = NA),
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
                 sourceURL = NA),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the parametrisation study area.",
                              "Note that `studyAreaLarge` is only used for parameter estimation, and",
                              "can be larger than the actual study area used for LandR simulations",
                              "(e.g., larger than `studyArea` in LandR Biomass_core)."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA),
    expectsInput("summaryPeriod", "numeric",
                 desc = "a numeric vector contains the start year and end year of summary",
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("vegTypeMap", "Raster", desc = NA)
  )
))

doEvent.LandWeb_output <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandWeb_output", "initialConditions",
                         eventPriority = 1)
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandWeb_output", "otherPlots",
                         eventPriority = 1)
    # sim <- scheduleEvent(sim, 0, "LandWeb_output", "allEvents", eventPriority = 7.5)
    sim <- scheduleEvent(sim, sim$summaryPeriod[1], "LandWeb_output", "allEvents",
                         eventPriority = 7.5)
  } else if (eventType == "initialConditions") {
    if (anyPlotting(P(sim)$.plots) && ("screen" %in% P(sim)$.plots)) {
      devCur <- dev.cur()
      ## if current plot dev is too small, open a new one
      if (is.null(dev.list())) {
        dev(x = devCur + 1, height = 7, width = 14)
        clearPlot()
      } else {
        if (dev.size()[2] < 14) {
          dev(x = devCur + 1, height = 7, width = 14)
          clearPlot()
        }
      }

      plotVTM(speciesStack = raster::mask(sim$speciesLayers, sim$studyAreaReporting) %>% raster::stack(),
              vegLeadingProportion = P(sim)$vegLeadingProportion,
              sppEquiv = sim$sppEquiv,
              sppEquivCol = P(sim)$sppEquivCol,
              colors = sim$sppColorVect,
              title = "Initial Types")

      dev(devCur)

      ## plot initial age map
      ageMap <- raster::mask(sim$standAgeMap, sim$studyAreaReporting) %>% raster::stack()
      Plot(ageMap, title = "Initial stand ages")
    }
  } else if (eventType == "allEvents") {
    if (time(sim) >= sim$summaryPeriod[1] && time(sim) <= sim$summaryPeriod[2]) {
      sim <- AllEvents(sim)
      sim <- scheduleEvent(sim,  time(sim) + P(sim)$summaryInterval,
                           "LandWeb_output", "allEvents", eventPriority = 7.5)
    }
  } else if (eventType == "otherPlots") {
    if (anyPlotting(P(sim)$.plots) && ("screen" %in% P(sim)$.plots)) {
      ## average age by FRI polygon
      mod$tsfOverTime <- ggPlotFn(sim$rstTimeSinceFire, sim$studyAreaReporting,
                                  sim$fireReturnInterval, current(sim)$eventTime, end(sim),
                                  mod$tsfOverTime, P(sim)$plotInitialTime, P(sim)$plotInterval,
                                  outputPath(sim))

      ## schedule future plots
      sim <- scheduleEvent(sim, times(sim)$current + P(sim)$.plotInterval, "LandWeb_output",
                           "otherPlots", eventPriority = 1)
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
  sim$vegTypeMap <- vegTypeMapGenerator(sim$cohortData, sim$pixelGroupMap,
                                        P(sim)$vegLeadingProportion,  mixedType = P(sim)$mixedType,
                                        sppEquiv = sim$sppEquiv, sppEquivCol = P(sim)$sppEquivCol,
                                        colors = sim$sppColorVect,
                                        doAssertion = getOption("LandR.assertions", TRUE))
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

  if (!suppliedElsewhere("summaryPeriod", sim))
    sim$summaryPeriod <- c(1000, 1500)

  if (!suppliedElsewhere("cohortData", sim))
    sim$cohortData <- data.table()

  if (!suppliedElsewhere("pixelGroupMap", sim))
    sim$pixelGroupMap <- raster()

  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath, cacheTags)
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    sim$sppEquiv <- LandR::sppEquivalencies_CA

    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## add default colors for species used in model
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol, newVals = "Mixed", palette = "Accent")
  }

  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    speciesLayersList <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyAreaLarge,
                               sppEquiv = sim$sppEquiv,
                               knnNamesCol = "KNN",
                               sppEquivCol = P(sim)$sppEquivCol,
                               # thresh = 10,
                               url = extractURL("speciesLayers"),
                               cachePath = cachePath(sim),
                               userTags = c(cacheTags, "speciesLayers"))
    #options(opts)

    writeRaster(speciesLayersList$speciesLayers,
                file.path(outputPath(sim), "speciesLayers.grd"),
                overwrite = TRUE)
    sim$speciesLayers <- speciesLayersList$speciesLayers
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
    sim$standAgeMap[] <- asInteger(sim$standAgeMap[])
  }

  return(invisible(sim))
}

ggPlotFn <- function(rstTimeSinceFire, studyAreaReporting, fireReturnInterval,
                     currTime, endTime, tsfOverTime, plotInitialTime, plotInterval, outPath) {
  tsfMap <- raster::mask(rstTimeSinceFire, studyAreaReporting)

  tsfDF <- data.table(tsf = tsfMap[], FRI = fireReturnInterval[]) %>% na.omit()
  tsfDF <- tsfDF[, list(
    time = as.numeric(currTime),
    meanAge = mean(tsf, na.rm = TRUE)), by = FRI]
  tsfDF[, FRI := factor(FRI)]

  tsfOverTime <- rbindlist(list(tsfOverTime, tsfDF))
  tsfOverTime <- tsfOverTime[!is.na(tsfOverTime$meanAge), ]

  if (length(unique(tsfOverTime$time)) > 1) {
    gg_tsfOverTime <- ggplot(tsfOverTime, aes(x = currTime, y = meanAge, col = FRI, ymin = 0)) +
      geom_line(size = 1.5) +
      theme(legend.text = element_text(size = 14))

    firstPlot <- isTRUE(currTime == plotInitialTime + plotInterval)
    title1 <- if (firstPlot) "Average age (TSF) by FRI polygon" else ""
    Plot(gg_tsfOverTime, title = title1, new = TRUE, addTo = "ageOverTime")

    if (currTime == endTime) {
      checkPath(file.path(outPath, "figures"), create = TRUE)
      ggsave(file.path(outPath, "figures", "average_age_(TSF)_by_FRI_polygon.png"), gg_tsfOverTime)
    }
  }
  return(tsfOverTime)
}
