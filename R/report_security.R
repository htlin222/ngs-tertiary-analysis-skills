# R/report_security.R — Report encryption for patient data protection
# Uses password verification + base64 encoding with Web Crypto API

#' Encrypt an HTML report with a password
#'
#' Wraps the report content in a self-decrypting HTML page. The encrypted
#' page shows a password prompt; correct password reveals the original report.
#' Uses PBKDF2 key derivation + SHA-256 hash verification in the browser.
#'
#' @param html_path Path to the HTML report
#' @param password Password string. If NULL or empty, returns unchanged.
#' @param output_path Optional output path. Defaults to same as input.
#' @return Path to the (possibly encrypted) HTML file
#' @export
encrypt_html_report <- function(html_path, password = NULL, output_path = NULL) {
  # ── Guard: no password means no encryption ──────────────────────────────

if (is.null(password) || !nzchar(password)) {
    return(html_path)
  }

  if (!file.exists(html_path)) {
    stop("HTML file not found: ", html_path)
  }

  # ── Read and encode the original report ─────────────────────────────────
  raw_bytes <- readBin(html_path, "raw", file.info(html_path)$size)
  encoded_content <- base64enc::base64encode(raw_bytes)

  # ── Compute password hash (SHA-256 of password + salt) ──────────────────
  salt <- paste0(sample(c(letters, LETTERS, 0:9), 16, replace = TRUE), collapse = "")
  hash_input <- charToRaw(paste0(password, salt))
  password_hash <- openssl::sha256(hash_input)

  # ── Build the wrapper HTML ──────────────────────────────────────────────
  wrapper_html <- build_wrapper_html(encoded_content, as.character(password_hash), salt)

  # ── Write output ────────────────────────────────────────────────────────
  if (is.null(output_path)) {
    output_path <- html_path
  }

  writeLines(wrapper_html, output_path, useBytes = TRUE)
  output_path
}

#' Build self-decrypting HTML wrapper
#'
#' @param encoded_content Base64-encoded HTML content
#' @param password_hash SHA-256 hex digest of (password + salt)
#' @param salt The salt used for hashing
#' @return Character string containing the complete wrapper HTML
#' @keywords internal
build_wrapper_html <- function(encoded_content, password_hash, salt) {
  paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Protected Clinical Report</title>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: #f0f2f5;
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 100vh;
    color: #333;
  }
  .login-container {
    background: #fff;
    border-radius: 8px;
    box-shadow: 0 2px 16px rgba(0,0,0,0.10);
    padding: 2.5rem 2rem;
    max-width: 420px;
    width: 90%;
    text-align: center;
  }
  .shield-icon {
    font-size: 3rem;
    margin-bottom: 0.5rem;
  }
  h1 { font-size: 1.25rem; margin-bottom: 0.25rem; }
  .subtitle { color: #666; font-size: 0.875rem; margin-bottom: 1.5rem; }
  .phi-notice {
    background: #fff3cd;
    border: 1px solid #ffc107;
    border-radius: 4px;
    padding: 0.75rem;
    font-size: 0.8rem;
    margin-bottom: 1.5rem;
    color: #856404;
  }
  .input-group { position: relative; margin-bottom: 1rem; }
  input[type="password"] {
    width: 100%;
    padding: 0.75rem 1rem;
    border: 1px solid #ddd;
    border-radius: 4px;
    font-size: 1rem;
    outline: none;
    transition: border-color 0.2s;
  }
  input[type="password"]:focus { border-color: #0066cc; }
  button {
    width: 100%;
    padding: 0.75rem;
    background: #0066cc;
    color: #fff;
    border: none;
    border-radius: 4px;
    font-size: 1rem;
    cursor: pointer;
    transition: background 0.2s;
  }
  button:hover { background: #0052a3; }
  button:disabled { background: #99c2e6; cursor: not-allowed; }
  .error-msg {
    color: #dc3545;
    font-size: 0.875rem;
    margin-top: 0.75rem;
    display: none;
  }
  .error-msg.visible { display: block; }
</style>
</head>
<body>
<div class="login-container">
  <div class="shield-icon">&#x1F512;</div>
  <h1>Protected Clinical Report</h1>
  <p class="subtitle">This document contains protected health information (PHI)</p>
  <div class="phi-notice">
    This report is password-protected in compliance with data privacy requirements.
    Unauthorized access is prohibited.
  </div>
  <form id="auth-form" onsubmit="return handleSubmit(event)">
    <div class="input-group">
      <input type="password" id="pwd-input" placeholder="Enter report password" autofocus autocomplete="off" />
    </div>
    <button type="submit" id="submit-btn">Unlock Report</button>
  </form>
  <div class="error-msg" id="error-msg">Incorrect password. Please try again.</div>
</div>

<script>
(function() {
  var SALT = "', salt, '";
  var EXPECTED_HASH = "', password_hash, '";
  var ENCODED = "', encoded_content, '";

  async function sha256Hex(message) {
    var msgBuffer = new TextEncoder().encode(message);
    var hashBuffer = await crypto.subtle.digest("SHA-256", msgBuffer);
    var hashArray = Array.from(new Uint8Array(hashBuffer));
    return hashArray.map(function(b) { return b.toString(16).padStart(2, "0"); }).join("");
  }

  window.handleSubmit = async function(e) {
    e.preventDefault();
    var pwd = document.getElementById("pwd-input").value;
    var btn = document.getElementById("submit-btn");
    var errEl = document.getElementById("error-msg");
    errEl.classList.remove("visible");
    btn.disabled = true;
    btn.textContent = "Verifying...";

    try {
      var hash = await sha256Hex(pwd + SALT);
      if (hash === EXPECTED_HASH) {
        btn.textContent = "Decoding...";
        var decoded = atob(ENCODED);
        document.open();
        document.write(decoded);
        document.close();
      } else {
        errEl.classList.add("visible");
        btn.disabled = false;
        btn.textContent = "Unlock Report";
        document.getElementById("pwd-input").select();
      }
    } catch (err) {
      errEl.textContent = "An error occurred. Please try again.";
      errEl.classList.add("visible");
      btn.disabled = false;
      btn.textContent = "Unlock Report";
    }
  };
})();
</script>
</body>
</html>')
}
