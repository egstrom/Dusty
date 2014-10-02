function prof_dusty, dir, scale, pobs, line, ps, res, w, prof_plot
  
; pobs the observed profile
; line the line number in the DUSTY image file (in the grid, 2 is
; 450mu and 3 is 850mu)
; ps the pixel scale
; res the spatial resolution (beam width)
; w the width of each annulus
; 
; This program will put the radial profile onto a 2D map, where the
; pixel size is ~beam/10. This will be done using interpolation. The
; map will then be convolved using a Gaussian (simulating the
; beam). Then the map is (once more) re-sampled to the same size as
; the observations, and a radial profile is extracted in exactly the
; same way as was done for the observations (using ANNULUS). Thus
; comparison is made easy.
;
;--- read in model results + initial exercises

  if file_test(dir+'mod_input.out') eq 1 then begin

     out = rascii(dir+'mod_input.out',skip=44,/double)
     r1 = out(3)*scale

     sc = float(ceil(10./(res/ps))) ; overall sampling factor - 
                                ; corresponds to ~beam/10
  
     robs = pobs(0,*)/ps        ; observations reported in ", change to pix
     srad = size(robs)

     pmod = fltarr(n_elements(pobs(0,*)))

     spec = rascii(dir+'mod_input.i001',skip=0,/double)
     tx = where(spec(line,*) le 10.D)
     tf = dblarr(n_elements(spec(*,0)),n_elements(tx))
     tf(*,*) = spec(*,tx)
     spec = tf

     b = spec(0,*)*r1
  
     rad = float(ceil(max(b)/ps)*sc) ; profile radius (in units of 
                                ; re-sampled pixels)
     dist_circle,dst,2*rad+1    ; following example by Jes, defining 
                                ; an array containing d from centre of array

     t = fltarr(n_elements(spec(line,*)))
     t(*) = spec(line,*)

     map = interpol([t,0.0],[b(*)*sc/ps,max(b(*)*sc/ps)*1.1],dst) 
                                ; interpolating to d

     map1 = filter_image(map,fwhm_gaussian=(res*sc/ps),/all_pixels)
     map2 = filter_image(map,fwhm_gaussian=(40.*sc/ps),/all_pixels)
     peak = 0.88*max(map1) + 0.12*max(map2)
     map = 0.88/peak*map1 + 0.12/peak*map2

     prof_plot = fltarr(2,rad+1)
     prof_plot(0,*) = findgen(rad+1)*ps/sc
     prof_plot(1,*) = map(rad,rad:2*rad)

     rad = rad/sc             ; re-re-sample onto original grid (units: pixels)
     map = congrid(map,2*rad+1,2*rad+1) ; using congrid

     for r=0,srad(2)-1 do begin ; doing the radial profile in the same 
                                ; manner as obs, also making sure obs 
                                ; do not extend beyond model results
        if robs(r)+w/2. le rad then pmod(r) = annulus(map,rad,rad,robs(r),w) $
        else pmod(r) = !values.f_nan
     endfor

                                ; added this little modification to
                                ; only look at model envelopes that
                                ; extend out to the same distance as
                                ; observed envelopes, i.e., adding a
                                ; criterion that the model envelope
                                ; must be at least as large as the
                                ; observed one.
;     if total(finite(pmod)) lt n_elements(pobs(0,*)) then pmod=!values.f_nan 

     return, pmod/pmod(0)    ; returning the normalized modelled radial profile

  endif else pmod=!values.f_nan

end
