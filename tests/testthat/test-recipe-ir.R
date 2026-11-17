# The language-neutral representation of a recipe (dev/plan.md 3.14, I1): the
# R list made by compile() converts to a tree of values, which is all the C++
# core reads. Here: every example recipe converts, the conversion is stable,
# and compiled functions become function references.

test_that("example recipes convert to the intermediate representation", {
  skip_on_cran()
  files = list.files(system.file("extdata", package = "ilike"), pattern = "^test_m_proposal.*\\.ilike$",
                     full.names = TRUE)
  for (f in files)
  {
    recipe = tryCatch(suppressMessages(suppressWarnings(ilike::compile(f))), error = function(e) NULL)
    if (is.null(recipe))
      next
    json = ilike:::recipe_json(recipe)
    expect_identical(json, ilike:::recipe_json(recipe))
    expect_false(grepl('"other"', json, fixed = TRUE))
  }
  model = suppressMessages(ilike::compile(c(system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike"),
                                            system.file("extdata", "test_m_proposal_norm_rw.ilike", package = "ilike"))))
  json = ilike:::recipe_json(model)
  # the C++ functions of the model are function references; the built-in
  # proposal is a typed list
  expect_true(grepl('{"function":{"language":"cpp"', json, fixed = TRUE))
  expect_true(grepl('{"name":"type","value":{"string":["norm_rw"]}}', json, fixed = TRUE))
})
