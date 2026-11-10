# Progress and status output. These go through message() so users can
# silence them with suppressMessages(); print methods still use cat().

# Like cat(..., sep = sep) but as a message (no newline is added).
bt_say <- function(..., sep = " ") {
  parts <- unlist(lapply(list(...), as.character))
  message(paste(parts, collapse = sep), appendLF = FALSE)
}

# Prints an object (e.g. a data frame) as a message.
bt_show <- function(x, ...) {
  message(paste(utils::capture.output(print(x, ...)), collapse = "\n"))
}
