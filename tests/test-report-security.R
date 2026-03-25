library(testthat)
library(here)

source(here("R/report_security.R"))

test_that("encrypt_html_report creates protected file", {
  tmp <- tempfile(fileext = ".html")
  writeLines("<html><body><h1>Test</h1><p>Patient data here</p></body></html>", tmp)
  output <- encrypt_html_report(tmp, password = "test123")
  expect_true(file.exists(output))
  content <- paste(readLines(output), collapse = "\n")
  # Should contain password prompt, not raw patient data
  expect_true(grepl("password", content, ignore.case = TRUE))
  expect_false(grepl("Patient data here", content))
  unlink(c(tmp, output))
})

test_that("encrypt_html_report with NULL password returns original", {
  tmp <- tempfile(fileext = ".html")
  writeLines("<html><body>test</body></html>", tmp)
  output <- encrypt_html_report(tmp, password = NULL)
  expect_equal(output, tmp)
  unlink(tmp)
})

test_that("encrypt_html_report with empty password returns original", {
  tmp <- tempfile(fileext = ".html")
  writeLines("<html><body>test</body></html>", tmp)
  output <- encrypt_html_report(tmp, password = "")
  expect_equal(output, tmp)
  unlink(tmp)
})

test_that("encrypt_html_report with custom output_path works", {
  tmp <- tempfile(fileext = ".html")
  out <- tempfile(fileext = "_protected.html")
  writeLines("<html><body>Secret data</body></html>", tmp)
  result <- encrypt_html_report(tmp, password = "pwd", output_path = out)
  expect_equal(result, out)
  expect_true(file.exists(out))
  # Original should still exist unchanged
  expect_true(file.exists(tmp))
  expect_true(grepl("Secret data", paste(readLines(tmp), collapse = "")))
  unlink(c(tmp, out))
})

test_that("encrypt_html_report errors on missing file", {
  expect_error(
    encrypt_html_report("/nonexistent/file.html", password = "test"),
    "not found"
  )
})

test_that("encrypted wrapper is self-contained HTML", {
  tmp <- tempfile(fileext = ".html")
  writeLines("<html><body>PHI content</body></html>", tmp)
  output <- encrypt_html_report(tmp, password = "secure")
  content <- paste(readLines(output), collapse = "\n")
  # Should be valid HTML with doctype
  expect_true(grepl("<!DOCTYPE html>", content))
  # Should contain the hash verification script
  expect_true(grepl("SHA-256", content))
  expect_true(grepl("EXPECTED_HASH", content))
  # Should contain PHI warning
  expect_true(grepl("protected health information", content, ignore.case = TRUE))
  # Should NOT contain the original content in plaintext

  expect_false(grepl("PHI content", content))
  unlink(c(tmp, output))
})
