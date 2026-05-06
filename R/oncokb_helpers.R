# R/oncokb_helpers.R — Convert HGVS.p 3-letter notation into the 1-letter
# short form OncoKB's /annotate/mutations endpoints actually parse.
#
# OncoKB's API matches alterations like "V600E" or "G12C". When given the
# 3-letter form ("Val600Glu" / "Gly12Ser"), the gene query is silently
# downgraded to "Other Biomarkers" and the response carries no oncogenicity
# or therapy data. TSO500's CombinedVariantOutput.tsv emits 3-letter
# notation by default — so every TSO500-driven OncoKB lookup must convert.

AA3_TO_AA1 <- c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
  Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
  Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
  Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
  Sec = "U", Pyl = "O", Asx = "B", Glx = "Z", Xle = "J",
  Xaa = "X", Ter = "*"
)

#' Convert a 3-letter HGVS.p alteration to OncoKB's expected 1-letter form.
#'
#' Examples:
#'   p.Gly12Ser           → G12S
#'   Val600Glu            → V600E
#'   p.(Lys17Arg)         → K17R
#'   Ala139CysfsTer89     → A139Cfs*89
#'   Tyr2148=             → Y2148=
#'   ?                    → (returned unchanged)
#'
#' Returns NA when input is NA or empty.
hgvsp_to_short <- function(hgvsp) {
  if (length(hgvsp) == 0) return(character())
  vapply(hgvsp, function(p) {
    if (is.na(p) || nchar(p) == 0) return(NA_character_)
    s <- p
    s <- sub("^p[.]", "", s)
    s <- gsub("[()]", "", s)

    # Replace every 3-letter AA token with its 1-letter form
    for (a3 in names(AA3_TO_AA1)) {
      s <- gsub(a3, AA3_TO_AA1[[a3]], s, fixed = TRUE)
    }
    s
  }, character(1), USE.NAMES = FALSE)
}
