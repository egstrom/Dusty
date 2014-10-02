function chi2_dusty, p850, p850_mod, p450, p450_mod, sed, sed_mod
  
; returning reduced chi2 for 850mu, 450mu and sed comparison
; separately (first three entries) and the combined chi2

  chi2sm = dblarr(5)

  if n_elements(p850_mod) le 1 then chi2sm(*)=!values.f_nan else begin
     if n_elements(p450_mod) le 1 then chi2sm(*)=!values.f_nan else begin
  
        np850 = n_elements(p850_mod)
        np450 = n_elements(p450_mod)
        nsed = n_elements(sed(0,*))

        for j=0,np450-1 do chi2sm(0) = chi2sm(0) + (p450(1,j)-p450_mod(j))^2/p450(2,j)^2
        for j=0,np850-1 do chi2sm(1) = chi2sm(1) + (p850(1,j)-p850_mod(j))^2/p850(2,j)^2
        for j=0,nsed-1 do chi2sm(2) = chi2sm(2) + total((sed(1,j)-sed_mod(j))^2/sed(2,j)^2)

        chi2sm(4) = total(chi2sm(0:1)) / (np850 + np450 - 3.)
        chi2sm(3) = total(chi2sm(0:2)) / (np850 + np450 + nsed - 3.)
        chi2sm(0) = chi2sm(0) / (np450 - 3.)
        chi2sm(1) = chi2sm(1) / (np850 - 3.)
        chi2sm(2) = chi2sm(2) / (nsed - 3.)
     endelse
  endelse

  return, chi2sm

end
