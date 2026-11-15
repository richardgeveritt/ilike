function my_dlnorm(tau,m,s,l)
  v = -log(tau[1]) - log(s) - 0.5*log(2*pi) - 0.5*((log(tau[1]) - m)/s)^2
  if (l==true)
    v
  else
    exp(v)
  end
end

function my_ones(n)
  ones(n)
end

function my_log_llhd(tau, y)
  sum(-0.5*log(2*pi) .- 0.5*(y .- tau[1]).^2)
end
