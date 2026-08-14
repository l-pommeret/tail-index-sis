dir.create("config", showWarnings = FALSE)

sensitivity <- expand.grid(
  model = c("A", "C", "D", "E"),
  signal = 1,
  n = 500L,
  p = 100L,
  rho = 0,
  reps = 30L,
  # Append the new value so legacy sensitivity cells retain config-row seeds.
  a = c(.25, .30, .35, .40, .45, .20),
  b = c(0, .10, .20, .30, .40),
  epsilon = .05,
  block_size = 32L,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
write.csv(sensitivity, "config/sensitivity.csv", row.names = FALSE)

scaling <- expand.grid(
  model = c("A", "C", "D", "E"), signal = 1,
  n = c(500L, 1000L, 2000L), p = 100L, rho = c(0, .5),
  reps = 100L, a = .20, b = .20, epsilon = .05, block_size = 64L,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
write.csv(scaling, "config/scaling.csv", row.names = FALSE)

main_dimension <- expand.grid(
  model="A", signal=1, n=c(500L,1000L,2000L),
  p=c(100L,500L,1000L,5000L), rho=0, reps=200L,
  a=.20, b=.20, epsilon=.05, block_size=128L,
  KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
write.csv(main_dimension, "config/main_dimension.csv", row.names=FALSE)

dependence <- expand.grid(
  model="A", signal=1, n=1000L, p=1000L,
  rho=c(0,.25,.5,.75), reps=200L, a=.20, b=.20, epsilon=.05,
  block_size=128L, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
write.csv(dependence, "config/dependence.csv", row.names=FALSE)

models <- rbind(
  data.frame(model="B", signal=c(1,.5,.25,.1)),
  data.frame(model=c("C","D","E","F","N"), signal=1))
models$n <- 1000L; models$p <- 1000L; models$rho <- 0; models$reps <- 200L
models$a <- .20; models$b <- .20; models$epsilon <- .05; models$block_size <- 128L
write.csv(models, "config/models.csv", row.names=FALSE)
