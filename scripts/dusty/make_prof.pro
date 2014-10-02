function make_prof, res, rad, w, p

  sc = 2.
  srad = n_elements(rad)
  b = findgen(240)*0.5+1
  r = 240
  pmod = fltarr(srad)

  dist_circle,dst,2*r+1         ; following example by Jes, defining 
                                  ; an array containing d from centre of array

  t = b^p

  map = interpol([t,0.0],[b(*),max(b(*))*1.1],dst) 
                                ; interpolating to d

  map1 = filter_image(map,fwhm_gaussian=(res*sc),/all_pixels)
  map2 = filter_image(map,fwhm_gaussian=(40.*sc),/all_pixels)
  map = 0.88*map1 + 0.12*map2
;  map = filter_image(map,fwhm_gaussian=(res*sc),/all_pixels)

  r = r/sc                ; re-re-sample onto original grid (units: pixels)
  map = congrid(map,2*r+1,2*r+1) ; using congrid

  for i=0,srad-1 do begin    ; doing the radial profile in the same 
                                ; manner as obs, also making sure obs 
                                ; do not extend beyond model results
     pmod(i) = annulus(map,r,r,rad(i),w)
  endfor

  pmod = pmod/pmod(0)
  prof=fltarr(3,srad)

  prof(0,*) = rad
  prof(1,*) = pmod
  prof(2,*) = 0.1*pmod

  return, prof               ; returning the normalized modelled radial profile

end
