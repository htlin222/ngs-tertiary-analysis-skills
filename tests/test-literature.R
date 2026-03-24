library(testthat)
library(here)
library(dplyr)

# Source the API clients
source(here("R/utils.R"))
source(here("R/api_clients.R"))

test_that("parse_pubmed_xml extracts articles", {
  # Create a minimal PubMed XML fixture
  xml_text <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2019//EN"
"https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID Version="1">12345678</PMID>
      <Article>
        <ArticleTitle>BRAF V600E mutations in melanoma: A comprehensive review</ArticleTitle>
        <Journal>
          <Title>Journal of Oncology</Title>
        </Journal>
        <PubDate>
          <Year>2023</Year>
        </PubDate>
        <AuthorList>
          <Author>
            <LastName>Smith</LastName>
            <Initials>JA</Initials>
          </Author>
          <Author>
            <LastName>Johnson</LastName>
            <Initials>MB</Initials>
          </Author>
        </AuthorList>
        <Abstract>
          <AbstractText>This study examines BRAF mutations.</AbstractText>
        </Abstract>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
  <PubmedArticle>
    <MedlineCitation>
      <PMID Version="1">87654321</PMID>
      <Article>
        <ArticleTitle>Targeted therapy in melanoma</ArticleTitle>
        <Journal>
          <Title>Cancer Research Today</Title>
        </Journal>
        <PubDate>
          <Year>2024</Year>
        </PubDate>
        <AuthorList>
          <Author>
            <LastName>Williams</LastName>
            <Initials>CD</Initials>
          </Author>
        </AuthorList>
        <Abstract>
          <AbstractText>A review of targeted therapies.</AbstractText>
        </Abstract>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>'

  result <- parse_pubmed_xml(xml_text)

  # Check structure
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)

  # Check columns exist
  expect_true("pmid" %in% names(result))
  expect_true("title" %in% names(result))
  expect_true("authors" %in% names(result))
  expect_true("journal" %in% names(result))
  expect_true("year" %in% names(result))
  expect_true("abstract" %in% names(result))

  # Check content
  expect_equal(result$pmid[1], "12345678")
  expect_equal(result$pmid[2], "87654321")
  expect_true(grepl("BRAF", result$title[1]))
  expect_true(grepl("melanoma", result$title[1], ignore.case = TRUE))
  expect_equal(result$year[1], "2023")
  expect_equal(result$year[2], "2024")

  # Check author parsing
  expect_true(grepl("Smith", result$authors[1]))
  expect_true(grepl("Johnson", result$authors[1]))
  expect_true(grepl("Williams", result$authors[2]))
})

test_that("parse_pubmed_xml handles empty XML", {
  xml_text <- '<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
</PubmedArticleSet>'

  result <- parse_pubmed_xml(xml_text)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("PubMed search tests skipped without key", {
  # Check if PUBMED_API_KEY is set
  api_key <- Sys.getenv("PUBMED_API_KEY", unset = "")

  skip_if(nchar(api_key) == 0, "PUBMED_API_KEY not set")
  skip_if(!requireNamespace("rentrez", quietly = TRUE), "rentrez package not installed")

  # If we reach here, we have required dependencies
  expect_true(nchar(api_key) > 0)
})
