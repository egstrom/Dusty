function chi2_sed, sed, sed_mod
  
; returning reduced chi2 for 850mu, 450mu and sed comparison
; separately (first three entries) and the combined chi2

  dum=0.D
  nsed = n_elements(sed(0,*))

  for j=0,nsed-1 do dum = dum + total((sed(1,j)-sed_mod(j))^2/sed(2,j)^2)

  chi2sm = dum / (nsed - 2.)

  return, chi2sm

end
