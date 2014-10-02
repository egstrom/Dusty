function mass_calcul, flux, dist, rad

  mass = 3.69e-6 * flux * dist^2 * (exp(16.94/10.)-1.) * (60./rad)^2

  return, mass

end
