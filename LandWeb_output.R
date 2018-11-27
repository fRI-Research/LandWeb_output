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
    defineParameter("summaryInterval", "numeric", 50, NA, NA, "This describes summary interval for this module"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "cohortData", objectClass = "data.table",
                 desc = "age cohort-biomass table hooked to pixel group map by pixelGroupIndex at
                 succession time step, this is imported from forest succession module",
                 sourceURL = NA),
    expectsInput("species", "data.table", "Columns: species, speciesCode, Indicating several features about species",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step, this is imported from
                 forest succession module",
                 sourceURL = NA),
    expectsInput(objectName = "summaryPeriod", objectClass = "numeric",
                 desc = "a numeric vector contains the start year and end year of summary",
                 sourceURL = NA),
    expectsInput(objectName = "vegLeadingProportion", objectClass = "numeric",
                 desc = "a number that define whether a species is lead for a given pixel",
                 sourceURL = NA)
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
  facVals <- pemisc::factorValues2(vtm, vtm[], att = "Species", na.rm = TRUE)
  df <- data.table(species = as.character(facVals), stringsAsFactors = FALSE)
  df <- df[!is.na(df$species)]
  df$species <- equivalentName(df$species, speciesEquivalency, "shortNames")
  df$cols <- equivalentName(df$species, speciesEquivalency, "cols")

  cols2 <- df$cols
  names(cols2) <- df$species
  initialLeadingPlot <- ggplot(data = df, aes(species, fill = species)) +
    scale_fill_manual(values = cols2) +
    geom_bar(position = "stack") +
    #labs(x = "Year", y = "Biomass by species") +
    theme(legend.text = element_text(size = 6), legend.title = element_blank(),
          axis.text = element_text(size = 6))

  Plot(initialLeadingPlot, title = c("Initial leading types"))
  Plot(vtm, title = "Initial leading types")
}

AllEvents <- function(sim) {
  sim$vegTypeMap <- vegTypeMapGenerator(sim$species, sim$cohortData, sim$pixelGroupMap,
                                        sim$vegLeadingProportion)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  if (!suppliedElsewhere("summaryPeriod", sim))
    sim$summaryPeriod <- c(1000, 1500)

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
      download.file(extractURL("species"), destfile = localSpeciesFilename)
    }
    sim$species <- read.csv(localSpeciesFilename, header = TRUE,
                            stringsAsFactors = FALSE) %>%
      data.table()
  }
  return(invisible(sim))
}
