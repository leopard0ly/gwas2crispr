## quiet Quarto invocation during CRAN checks
## This .onLoad hook prevents accidental invocation of Quarto during package loading when
## the environment variable `_GWAS2CRISPR_DISABLE_QUARTO` is set. On CRAN,
## tests set this variable to "1" so that any call to system2("quarto", ...) does
## nothing. The original system2 is restored otherwise.

.onLoad <- function(libname, pkgname) {
  # only override system2 if a special environment variable is set
  disable_quarto <- Sys.getenv("_GWAS2CRISPR_DISABLE_QUARTO", unset = NA_character_)
  if (!is.na(disable_quarto) && nzchar(disable_quarto)) {
    # capture original system2
    orig_system2 <- base::system2
    # unlock binding to allow assignment
    unlockBinding("system2", asNamespace("base"))
    base::system2 <- function(command, args = character(), ...) {
      # if the command is 'quarto', skip execution and return success code
      if (identical(command, "quarto")) {
        return(invisible(0L))
      }
      orig_system2(command, args, ...)
    }
    # relock binding after override
    lockBinding("system2", asNamespace("base"))
  }
}
