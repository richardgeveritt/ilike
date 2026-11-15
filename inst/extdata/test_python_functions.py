import numpy as np

def my_ones(n):
  return np.ones(n)

def my_log_prior(tau):
  t = float(np.asarray(tau).reshape(-1)[0])
  return -np.log(t) - 0.5*np.log(2*np.pi) - 0.5*np.log(t)**2

def my_log_llhd(parameters, data):
  y = np.asarray(data["y"]).reshape(-1)
  tau = float(np.asarray(parameters["tau"]).reshape(-1)[0])
  return float(np.sum(-0.5*np.log(2*np.pi) - 0.5*(y - tau)**2))
