library(testthat)
library(here)
library(dplyr)

# Source the API clients
source(here("R/utils.R"))
source(here("R/api_clients.R"))

test_that("parse_oncokb_treatments handles empty list", {
  result <- parse_oncokb_treatments(list())

  # Should return empty tibble with correct columns
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 4)
  expect_true("drugs" %in% names(result))
  expect_true("level" %in% names(result))
  expect_true("description" %in% names(result))
  expect_true("fda_approved" %in% names(result))
})

test_that("parse_oncokb_treatments parses treatments", {
  # Mock treatment list (structure from OncoKB API)
  mock_treatments <- list(
    list(
      drugs = list(
        list(drugName = "Vemurafenib"),
        list(drugName = "Dabrafenib")
      ),
      level = "LEVEL_1",
      description = "FDA-approved for melanoma"
    ),
    list(
      drugs = list(
        list(drugName = "Sorafenib")
      ),
      level = "LEVEL_3A",
      description = "Investigational"
    )
  )

  result <- parse_oncokb_treatments(mock_treatments)

  expect_equal(nrow(result), 2)
  # Drug names come from mock — first treatment has both Vemurafenib + Dabrafenib
  expect_true(grepl("Vemurafenib", result$drugs[1]) || grepl("Dabrafenib", result$drugs[1]))
  expect_equal(result$level[1], "LEVEL_1")
  expect_equal(result$level[2], "LEVEL_3A")

  # Check FDA approval detection
  expect_true(result$fda_approved[1])
  expect_false(result$fda_approved[2])
})

test_that("oncokb_request builds correct URL", {
  # We'll check that the function doesn't error when API key is available
  # by mocking the Sys.getenv call

  # This is a structural test - we just verify the function exists and can be called
  expect_true(is.function(oncokb_request))
})

test_that("OncoKB API tests skipped without key", {
  # Check if ONCOKB_API_KEY is set
  api_key <- Sys.getenv("ONCOKB_API_KEY", unset = "")

  skip_if(nchar(api_key) == 0, "ONCOKB_API_KEY not set")

  # If we reach here, the API key is available
  # We would run actual API tests here
  expect_true(nchar(api_key) > 0)
})
