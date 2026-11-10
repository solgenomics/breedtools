# Runs `expr`, silences its messages and returns them as lines of text,
# like capture.output() does for printed output.
capture_msgs <- function(expr) {
  msgs <- character()
  withCallingHandlers(expr, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  out <- unlist(strsplit(paste(msgs, collapse = ""), "\n", fixed = TRUE))
  if (is.null(out)) character() else out
}

# agricolae's design functions switch the RNG kind to "Super-Duper" and leave
# it there, which changes the random data of later tests. Restores the kind
# and the seed when the calling test ends.
local_agricolae_rng <- function(env = parent.frame()) {
  kind <- RNGkind()
  withr::local_preserve_seed(.local_envir = env)
  withr::defer(do.call(RNGkind, as.list(kind)), envir = env)
}
