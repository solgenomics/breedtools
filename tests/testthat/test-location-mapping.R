test_that("read_locations standardizes coordinate files", {
  locations <- data.frame(
    location_name = c("Loc1", "Loc2"),
    lat = c("35.1", "36.2"),
    lon = c("-78.9", "-79.4")
  )

  result <- read_locations(
    locations,
    reference_locations = c("loc1", "loc2"),
    ignore_case = TRUE
  )

  expect_equal(names(result), c("location", "latitude", "longitude"))
  expect_equal(result$location, c("Loc1", "Loc2"))
  expect_type(result$latitude, "double")
  expect_type(result$longitude, "double")
})

test_that("read_locations requires matching reference locations", {
  locations <- data.frame(
    location = "Loc1",
    latitude = 35.1,
    longitude = -78.9
  )

  expect_error(
    read_locations(locations, reference_locations = c("Loc1", "Loc2")),
    "Missing coordinates"
  )
})

test_that("read_locations validates against phenotype data automatically", {
  locations <- data.frame(
    location = c("loc1", "loc2"),
    latitude = c(35.1, 36.2),
    longitude = c(-78.9, -79.4)
  )
  pheno <- data.frame(
    locationName = c("Loc1", "Loc2", "Loc1"),
    germplasmName = c("G1", "G1", "G2"),
    yield = c(10, 12, 11)
  )

  result <- read_locations(
    locations,
    phenotype_data = pheno,
    ignore_case = TRUE
  )

  expect_equal(nrow(result), 2)
})

test_that("read_locations validates against stratification result fields", {
  locations <- data.frame(
    location = c("Env1", "Env2"),
    latitude = c(35.1, 36.2),
    longitude = c(-78.9, -79.4)
  )
  stratification_result <- list(
    env_info = data.frame(
      environment = c("Env1", "Env2"),
      location = c("1", "2"),
      environment_label = c("Location One", "Location Two")
    )
  )

  result <- read_locations(
    locations,
    phenotype_data = stratification_result,
    match_by = "environment"
  )

  expect_equal(result$location, c("Env1", "Env2"))
})

test_that("map_environment_groups returns a leaflet map", {
  stratification_result <- list(
    group_membership = data.frame(
      group_id = c("Group_1", "Group_1"),
      environment = c("Env1", "Env2"),
      location = c("Loc1", "Loc2"),
      environment_label = c("Loc1 / Trial1", "Loc2 / Trial1")
    ),
    ungrouped = data.frame(
      environment = "Env3",
      location = "Loc3",
      environment_label = "Loc3 / Trial1"
    )
  )

  locations <- read_locations(data.frame(
    location = c("loc1", "loc2", "loc3"),
    latitude = c(35.1, 36.2, 37.3),
    longitude = c(-78.9, -79.4, -80.1)
  ))

  result <- map_environment_groups(
    stratification_result,
    locations,
    ignore_case = TRUE
  )

  expect_s3_class(result, "leaflet")
  expect_s3_class(result, "htmlwidget")
})

test_that("map_environment_groups can match coordinates by environment", {
  stratification_result <- list(
    group_membership = data.frame(
      group_id = c("Group_1", "Group_1"),
      environment = c("Env1", "Env2"),
      location = c("1", "2"),
      environment_label = c("Location One", "Location Two")
    ),
    ungrouped = data.frame(
      environment = character(),
      location = character(),
      environment_label = character()
    )
  )

  locations <- read_locations(data.frame(
    location = c("env1", "env2"),
    latitude = c(35.1, 36.2),
    longitude = c(-78.9, -79.4)
  ))

  result <- map_environment_groups(
    stratification_result,
    locations,
    ignore_case = TRUE,
    match_by = "environment"
  )

  expect_s3_class(result, "leaflet")
})

test_that("map_environment_groups uses open, keyless maps by default", {
  pheno_file <- system.file("extdata", "breedbase_pheno.csv", package = "breedtools")
  loc_file <- system.file("extdata", "breedbase_locations.csv", package = "breedtools")
  result <- run_environment_stratification(pheno_file, trait = "fresh root yield")
  locations <- read_locations(loc_file, phenotype_data = result)

  map <- map_environment_groups(result, locations)
  methods <- vapply(map$x$calls, `[[`, character(1), "method")
  tile_urls <- vapply(map$x$calls[methods == "addTiles"], function(x) x$args[[1]], character(1))
  expect_length(tile_urls, 3)
  expect_true(all(grepl("openstreetmap\\.org|opentopomap\\.org|maps\\.eox\\.at", tile_urls)))
  expect_false(any(grepl("arcgis|esri|google|mapbox", tile_urls, ignore.case = TRUE)))
  expect_false(any(grepl("key|token", tile_urls, ignore.case = TRUE)))
  control <- map$x$calls[[which(methods == "addLayersControl")]]
  expect_equal(control$args[[1]], c("Street", "Topographical", "Satellite"))

  sat <- map_environment_groups(result, locations, tiles = "satellite")
  expect_match(sat$x$calls[[1]]$args[[1]], "s2cloudless")
  old_name <- map_environment_groups(result, locations, tiles = "breedbase")
  expect_identical(old_name$x$calls[[1]]$args[[1]], tile_urls[[1]])
  expect_error(map_environment_groups(result, locations, tiles = ""), "one map name")
})
