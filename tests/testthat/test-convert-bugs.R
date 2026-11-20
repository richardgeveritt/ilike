# Converting BUGS, JAGS and NIMBLE models to the ilike language (dev/plan.md
# 3.4, C231): the translated models' log densities against R's, and refusals.

bugs_value = function(model, data, P, dialect = "jags")
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  ilike::convert_bugs_model(model, data, dialect = dialect, output = f)
  ilike::evaluate_log(ilike::compile(f), parameters = P)
}

test_that("JAGS models with precisions, link functions, inprod and truncation", {
  y = c(1.2, 0.4, 2.1, 1.7); x = c(-1, 0, 1, 2); z = c(1, 0, 1, 1)
  model = "model {
    for (i in 1:N) { y[i] ~ dnorm(mu[i], tau)   # likelihood
                     mu[i] <- alpha + beta * x[i] }
    alpha ~ dnorm(0, 0.001)
    beta ~ dnorm(0, 1.0E-2)
    tau ~ dgamma(0.1, 0.1)
  }"
  text = ilike::convert_bugs_model(model, list(N = 4, y = y, x = x))
  expect_match(text, "move m = nuts([alpha, beta, tau]", fixed = TRUE)
  expect_equal(bugs_value(model, list(N = 4, y = y, x = x), list(alpha = 0.3, beta = 0.6, tau = 1.7)),
               dnorm(0.3, 0, 1/sqrt(0.001), log = TRUE) + dnorm(0.6, 0, 10, log = TRUE) + dgamma(1.7, 0.1, 0.1, log = TRUE) +
                 sum(dnorm(y, 0.3 + 0.6*x, 1/sqrt(1.7), log = TRUE)), tolerance = 1e-10)
  X = cbind(x, x^2)
  expect_equal(bugs_value("model {
      for (i in 1:N) { logit(p[i]) <- a + inprod(X[i, ], b[]); z[i] ~ dbern(p[i]) }
      a ~ dnorm(0, 0.25) T(0, )
      for (k in 1:2) { b[k] ~ dt(0, 1, 3) }
    }", list(N = 4, z = z, X = X), list(a = 0.5, b = c(0.2, -0.3))),
    dnorm(0.5, 0, 2, log = TRUE) - pnorm(0, 0, 2, lower.tail = FALSE, log.p = TRUE) + sum(dt(c(0.2, -0.3), 3, log = TRUE)) +
      sum(dbinom(z, 1, plogis(0.5 + X %*% c(0.2, -0.3)), log = TRUE)), tolerance = 1e-10)
})

test_that("NIMBLE's named parameterisations, and an R expression as the model", {
  y = c(1.2, 0.4, 2.1, 1.7)
  expect_equal(bugs_value(quote({
      mu ~ dnorm(0, sd = 10)
      sigma ~ dunif(0, 5)
      for (i in 1:N) y[i] ~ dnorm(mu, sd = sigma)
      w ~ dweib(2, scale = 3)
      v ~ dgamma(2, scale = 0.5)
    }), list(N = 4, y = y, w = 1.3, v = 0.7), list(mu = 1, sigma = 0.8), dialect = "nimble"),
    dnorm(1, 0, 10, log = TRUE) + dunif(0.8, 0, 5, log = TRUE) + sum(dnorm(y, 1, 0.8, log = TRUE)) + dweibull(1.3, 2, 3, log = TRUE) +
      dgamma(0.7, 2, scale = 0.5, log = TRUE), tolerance = 1e-10)
})

test_that("vector distributions: dmnorm with a precision matrix as data, dmulti, ddirch, and JAGS's dweib", {
  Om = matrix(c(2, 0.5, 0.5, 1), 2)
  model = "model { m[1:2] ~ dmnorm(m0[], Omega[,]); t ~ dweib(2, 0.5); c[1:3] ~ dmulti(p[], 6); p[1:3] ~ ddirch(a[]) }"
  data = list(m0 = c(0, 0), Omega = Om, t = 1.1, c = c(1, 3, 2), a = c(1, 2, 1))
  text = ilike::convert_bugs_model(model, data)
  expect_match(text, "Omega_inverse", fixed = TRUE)
  expect_match(text, "move m = nuts([m, p]", fixed = TRUE)
  S = solve(Om); m = c(0.3, -0.2)
  expect_equal(bugs_value(model, data, list(m = m, p = c(0.2, 0.5, 0.3))),
               -log(2*pi) - 0.5*determinant(S)$modulus[1] - 0.5*sum(m*solve(S, m)) + dweibull(1.1, 2, 0.5^(-1/2), log = TRUE) +
                 dmultinom(c(1, 3, 2), 6, c(0.2, 0.5, 0.3), log = TRUE) + lgamma(4) - lgamma(2) + log(0.5), tolerance = 1e-10)
})

test_that("what cannot be translated is refused with a message", {
  refuse = function(model, data, message) expect_error(ilike::convert_bugs_model(model, data), message, fixed = TRUE)
  refuse("model { for (i in 1:2) { y[i] ~ dnorm(mu, 1) }; mu ~ dnorm(0, 1) }", list(y = c(1, NA)), "has missing values")
  refuse("model { y ~ dinterval(mu, 1); mu ~ dnorm(0, 1) }", list(y = 1), "censoring (dinterval)")
  refuse("model { m[1:2] ~ dmnorm(m0[], Om[,]); Om[1:2, 1:2] ~ dwish(R[,], 3) }", list(m = c(1, 2), m0 = c(0, 0), R = diag(2)),
         "must be data")
  refuse("model { probit(p) <- a; a ~ dnorm(0, 1); y ~ dbern(p) }", list(y = 1), "the link function probit")
  refuse("model { y ~ dnorm(mu, 1); mu ~ dgamma(1, 1) T(0, 2) }", list(y = 1), "truncation T(,) is translated for dnorm only")
  refuse("model { y ~ dnorm(step(mu), 1); mu ~ dnorm(0, 1) }", list(y = 1), "the function step cannot be translated")
})
