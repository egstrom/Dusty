pro env_dusty, dir, home, name, alpham, Ym, taum, scale, D

  out = rascii(dir+'mod_input.out',skip=44,/double)
  t = rascii(dir+'mod_input.rtb',/double)

  r1 = out(3)*scale

;--- calculating envelope properties

; different constants
  k100 = 86.5 ; Dust opacity from Ossenkopf & Henning (1994) [cm2/g]
  mu = 2.80                     ; mean molecular mass
  mH = 1.6733e-24               ; H mass
  fr = 0.01                     ; dust/gas
  MS = 1.989e33                 ; Solar mass
  au = 1.496e13                 ; AU in cm

; first transforming best-fit results into numbers

  p = fltarr(1)
  p(*) = alpham
  Y = fltarr(1)
  Y(*) = Ym
  tau = fltarr(1)
  tau(*) = taum
  ri = r1*D*au

  y10 = interpol(t(0,*),t(5,*),10.)
  if y10 gt Y then y10=Y

  NH2 = tau /(k100*fr*mu*mH) ; Eqns. 1+2 from Schoier et al. (2002)

  if p ne 1. then begin
     ni = NH2/ri*(1.-p)/(Y^(1.-p)-1.)
     NH2_10 = ni*ri/(1.-p)*(y10^(1.-p)-1.)
     MH2 = 4.*!pi*mu*mH*ni*double(ri)^3/(3.-p)*(y10^(3.-p)-1.)/MS
  endif else begin
     ni = NH2/ri/alog(Y)
     MH2 = 4.*!pi*mu*mH*ni*double(ri)^3/(3.-p)*(y10^(3.-p)-1.)/MS
     NH2_10 = ni*ri*alog(y10)
  endelse

  x = where((t(5,*)-10.) eq min(t(5,*)-10.))

  openw,unit,home+'/sources/'+name+'/res_'+name+'.dat',/get_lun

  printf,unit,'----------------------------------------'
  printf,unit,'Best fit model results are:'
  printf,unit,['                r1 [AU] = ',string(format='(1X,F7.2)',ri/au)]
  printf,unit,['            r(10K) [AU] = ',string(format='(e8.2)',ri*y10/au)]
  printf,unit,['              ni [cm-3] = ',string(format='(e8.2)',ni)]
  printf,unit,['      n(1000 AU) [cm-3] = ',string(format='(e8.2)',ni*(1.496e16/ri)^(-p))]
  printf,unit,['         n(10 K) [cm-3] = ',string(format='(e8.2)',ni*y10^(-p))]
  printf,unit,['     N(H2)(10 K) [cm-2] = ',string(format='(e8.2)',NH2_10)]
  printf,unit,['Envelope mass(10K) [MS] = ',string(format='(1X,F7.3)',MH2)]

  close, unit
  free_lun, unit

end
