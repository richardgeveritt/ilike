# Rook neighbour lists on an L x L lattice, for the tests of models with
# neighbour sums (test-structured.R) and unnormalised distributions
# (test-unnormalised.R)
lattice = function(L)
{
  # rook neighbours on an L x L lattice, padded with 0
  site = function(r, c) (c - 1)*L + r
  nb_list = lapply(seq_len(L*L), function(s) {
    r = (s - 1) %% L + 1; c = (s - 1) %/% L + 1
    out = c()
    if (r > 1) out = c(out, site(r - 1, c)); if (r < L) out = c(out, site(r + 1, c))
    if (c > 1) out = c(out, site(r, c - 1)); if (c < L) out = c(out, site(r, c + 1))
    out })
  nn = sapply(nb_list, length)
  nb = t(sapply(nb_list, function(v) c(v, rep(0, 4 - length(v)))))
  list(n = L*L, nb_list = nb_list, nn = nn, nb = nb,
       data = paste0("nn = [", paste(nn, collapse = ", "), "]; nb = [", paste0("[", apply(nb, 1, paste, collapse = ", "), "]", collapse = ", "), "]"))
}
