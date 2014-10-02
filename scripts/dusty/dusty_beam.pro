function dusty_beam, dir, sed, beam, scale, y

  out = rascii(dir+'mod_input.out',skip=44,/double)
  r1 = out(3)*scale

  wave = [25.00, 50.00, 60.00, 70.00, 100.00, 160.00, 350.00, 450.00, 750.00, 800.00, 850.00, 1100.00, 1300.00, 2000.00, 3000.00]

  ps = 0.25
  res = beam

  rad = ceil(float(y)*r1/ps)
  dist_circle,dst,2*rad+1

  spec = rascii(dir+'mod_input.i001',skip=0,/double)
  result = fltarr(n_elements(beam))

  for i=0,n_elements(beam)-1 do begin
     line = where(wave eq sed(0,i))+2
     tx = where(spec(line,*) le 10.D)
     tf = dblarr(n_elements(spec(*,0)),n_elements(tx))
     tf(*,*) = spec(*,tx)
     spec = tf

     t = fltarr(n_elements(spec(line,*)))
     t(*) = spec(line,*)

     b = spec(0,*)*r1/ps

     map = interpol([t,0.0],[b(*),max(b(*))*1.1],dst) 
                                ; interpolating to d

     psf = exp(-2.77*dst^2/(beam(i)/ps)^2)
     map = map*psf

     result(i) = total(map)*ps^2*1e-30*(3e10/wave(line-2)*1e4)
  endfor

  return,result

end
