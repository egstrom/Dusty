pro f850_extr, rad

  home = '/data/kristensen/wish/dusty_model/'

  data = rascii(home+'scuba_position.txt', header = in_name)
  temp = rascii(home+'name.txt', header = out_name)
  dist = rascii(home+'distance.txt')

  stot = n_elements(in_name)
;  rad = 3.33                    ; radius in pixels to be used
  ps = 6.                       ; pixel-scale

  flux = fltarr(stot)
  mass = fltarr(stot)

  for i=0,stot-1 do begin
     if data(0,i) ne 0. then begin
        in=in_name(i)+'_850um'
        flux(i)=flux_extr(in, data(0,i), data(1,i), 2, rad/2., data(2:5,i))

        mass(i) = mass_calcul(flux(i), dist(i), rad*ps)
     endif else begin
        flux(i) = 0
        mass(i) = 0
     endelse
  endfor

;  flux = flux * !pi * rad^2.

  for i=0,stot-1 do begin
     print,out_name(i), flux(i), mass(i)
  endfor

end
