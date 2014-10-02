function sed_dusty, dir, sed, intf

; Interpolates modelled SED to the wavelengths of the observed
; SED. Returns an array containing the modelled SED for all values of
; tau.

  sed_mod_t = rascii(dir+'mod_input.stb',/double)

  sed_mod = dblarr(n_elements(sed(0,*)))

  sed_mod(*) = hermite(sed_mod_t(0,*),sed_mod_t(1,*),sed(0,*))
                                ; hermite-interpolation (astro-lib)
                                ; proved to give much better results
                                ; than the internal interpol routine
                                ; in idl independent of the choice of
                                ; interpolation (linear, lsquadratic,
                                ; quadratic, spline)

  sed_mod = sed_mod * intf

  for i=0,n_elements(sed_mod(*))-1 do begin
     if finite(sed_mod(i)) ne 1 then print,'SED problems!'
  endfor

  return, sed_mod

end
