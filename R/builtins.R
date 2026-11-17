# The registry of built-ins (dev/plan.md 3.14, I2), as a data frame with one
# row per argument of each built-in: the source for the documentation's
# argument lists, and used by the tests. Internal (ilike:::builtin_table()).
builtin_table <- function()
{
  registry = builtin_registry_table()
  rows = list()
  for (entry in registry)
  {
    roles = paste(entry$roles, collapse = ", ")
    if (length(entry$arguments) == 0)
    {
      rows[[length(rows) + 1]] = data.frame(builtin = entry$name, roles = roles,
                                            argument = NA_character_, type = NA_character_,
                                            required = NA, default = NA_character_,
                                            description = entry$description,
                                            stringsAsFactors = FALSE)
      next
    }
    for (argument in entry$arguments)
    {
      rows[[length(rows) + 1]] = data.frame(builtin = entry$name, roles = roles,
                                            argument = argument$name, type = argument$type,
                                            required = argument$required, default = argument$default,
                                            description = argument$description,
                                            stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
