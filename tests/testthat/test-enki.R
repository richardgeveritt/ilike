# Iterative ensemble Kalman inversion (IENKI) on a conjugate model:
# theta ~ N(0,1), y = 1 ~ N(theta,1). At power phi the tempered posterior is
# N(phi/(1+phi), 1/(1+phi)), and log Z(phi) = -(phi/2) log(2 pi)
# - log(1+phi)/2 - phi/(2(1+phi)); the posterior is N(0.5, 0.5).

enki_log_z = function(phi) -(phi/2)*log(2*pi) - 0.5*log(1+phi) - phi/(2*(1+phi))

enki_run = function(sequence_lines = character(0), members = 4000)
{
  f = tempfile(fileext = ".ilike")
  writeLines(c("/***prior,ilike::norm(theta,0,1)***/", "",
               "/***data***/", "Data set_data() { return Data(\"y\", 1.0); }", "",
               "/***linear_gaussian_data_model,ilike::linear_gaussian_data_model(y,theta,1,1)***/", "",
               sequence_lines), f)
  name = basename(tempfile("enki_"))
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  z = ilike::IENKI(ilike::compile(f), results_name = name, number_of_ensemble_members = members, seed = 1)
  list(z = z, output = ilike::load_enk_output(name))
}

test_that("IENKI moves the ensemble to the posterior and estimates the evidence", {
  skip_on_cran()
  # one step, from the prior to the posterior (the default sequence)
  run = enki_run()
  expect_equal(run$z, enki_log_z(1), tolerance = 0.02)
  final = run$output[run$output$Target == max(run$output$Target), ]
  expect_equal(mean(final$Value), 0.5, tolerance = 0.1)
  expect_equal(var(final$Value), 0.5, tolerance = 0.1)

  # three steps (before the fix of C182 the ensemble never moved, so every
  # stored ensemble was the prior and the later increments were wrong)
  run = enki_run("/***enk_sequence,ilike::annealing(power,c(0,0.25,0.5,1))***/")
  expect_equal(run$z, enki_log_z(1), tolerance = 0.02)
  for (phi in c(0.25, 0.5, 1))
  {
    step = run$output[run$output$TargetParameters == paste0("power=", phi, ";"), ]
    expect_equal(mean(step$Value), phi/(1+phi), tolerance = 0.1)
    expect_equal(var(step$Value), 1/(1+phi), tolerance = 0.1)
    expect_equal(step$NormalisingConstant[1], enki_log_z(phi), tolerance = 0.02)
  }
})
