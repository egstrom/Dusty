function chi2_profile, p850, p850_mod, p450, p450_mod
  
; returning reduced chi2 for 850mu, 450mu and sed comparison
; separately (first three entries) and the combined chi2

  dum = dblarr(2)

  if n_elements(p850_mod) le 1 then chi2sm=!values.f_nan else begin
     if n_elements(p450_mod) le 1 then chi2sm=!values.f_nan else begin
  
        np850 = n_elements(p850_mod)
        np450 = n_elements(p450_mod)

        for j=0,np450-1 do dum(0) = dum(0) + (p450(1,j)-p450_mod(j))^2/p450(2,j)^2
        for j=0,np850-1 do dum(1) = dum(1) + (p850(1,j)-p850_mod(j))^2/p850(2,j)^2

        chi2sm = total(dum(0:1)) / (np850 + np450 - 1.)
     endelse
  endelse

  return, chi2sm

end
