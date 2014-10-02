pro result_dusty_beam, name, dist, lum, chi2

; This is the top-level program that controls the comparison of dust
; continuum observations with DUSTY.
;
; Radial profiles should be given in arrays, where the first column
; contains the radial coordinate ["], second column the
; normalized brightness (normalized to the first point) and the third
; column the error
;
; The SED should be given in an array, where the first column contains
; the wavelength [mu], the second the flux [Jy] and the third the
; error.
;
; If the directory does not exist, the function returns chi2 = NaN


  tid = systime(1)              ; just a small timer

;--- read in observations

  home = '/data/kristensen/wish/dusty_model/'

  sed_t = rascii(home+'sed/sed_'+name+'.txt',/double)
  st = size(sed_t)

  sed = dblarr(3,st(2))
  sed(0:1,*) = sed_t(0:1,*)
  dum1 = where(sed(0,*) eq 450.)
  dum2 = where(sed(0,*) eq 850.)
  sed(2,*)=sed(1,*)*0.3         ; assume 20% error on all sed values
                                ; following Jorgensen et al. (2002)
  sed(2,dum1)=sed(1,dum1)*0.05
  sed(2,dum2)=sed(1,dum2)*0.01
  sed(1,*)=sed(1,*)*1e-30*(3e10/sed(0,*)*1e4) ; change units from Jy to Wcm-2
  sed(2,*)=sed(2,*)*1e-30*(3e10/sed(0,*)*1e4)

  bm = where(sed_t(2,*) ne 0.)
  nobm = where(sed_t(2,*) eq 0.)
  beam = sed_t(2,bm)

;--- object specific parameters/constants

  D = dist                      ; absolute distance to object [pc]
  Lb = lum                      ; absolute bolometric luminosity [solar units]
  scale = sqrt(Lb/1e4)/1.496e13/D ; distance scale factor (scales as L^0.5)
  intf = Lb*double(3.839e26) / ((D*double(3.086e18))^2*4D*!pi)
                                ; integrated bolometric flux [W cm-2]

;--- Instrument specific parameters

  res850 = 19.5                 ; spatial resolution at 850 mu ["]
  ps850 = 6.                    ; pixel scale ["/pix]
  
  res450 = 11.0
  ps450 = 3.

  w850 = 1.5 ; distance between two points on the observed radial profile [pix]
  w450 = 2.0

  l850 = 12                      ; column to be read in DUSTY image file
  l450 = 9

  p = -1.5                      ; power-law exp to be used if no profile exists

  if file_test(home+'profiles/prof_'+name+'_850mu.txt') eq 1 then begin
     p850 = rascii(home+'profiles/prof_'+name+'_850mu.txt',/double) 
     p850(0,*) = ps850*p850(0,*)
     endif else $
        p850 = make_prof(res850, [1.5,3,4.5,6,7.5,9]*ps850, w850*ps850, p)

  if file_test(home+'profiles/prof_'+name+'_450mu.txt') eq 1 then begin
     p450 = rascii(home+'profiles/prof_'+name+'_450mu.txt',/double) 
     p450(0,*) = ps450*p450(0,*)
     endif else $
        p450 = make_prof(res450, [2,4,6,8,10,12,14,16,18]*ps450, w450*ps450, p)

;--- define model grid (to be used when constructing directory names below)

; refers specifically to the grid in LEKs home directory

  dusty_dir = '/data1/kristensen/home/scuba/dusty/results/'

  alpha=['0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4',$
         '1.5','1.6','1.7','1.8','1.9','2.0','2.1','2.2','2.3']
  Y = ['500','600','700','800','900','1000','1100','1200','1300','1400',$
       '1500','1600','1700','1800','1900','2000','2100','2200','2300',$
       '2400','2500','2600','2700','2800','2900','3000']
  tau = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1',$
         '1.2','1.3','1.4','1.5','1.6','1.7','1.8','1.9','2.0','2.1',$
         '2.2','2.3','2.4','2.5','2.6','2.7','2.8','2.9','3.0','3.1',$
         '3.2','3.3','3.4','3.5','3.6','3.7','3.8','3.9','4.0','4.1',$
         '4.2','4.3','4.4','4.5','4.6','4.7','4.8','4.9','5.0','5.1',$
         '5.2','5.3','5.4','5.5','5.6','5.7','5.8','5.9','6.0','6.1',$
         '6.2','6.3','6.4','6.5','6.6','6.7','6.8','6.9','7.0','7.1',$
         '7.2','7.3','7.4','7.5','7.6','7.7','7.8','7.9','8.0','8.1',$
         '8.2','8.3','8.4','8.5','8.6','8.7','8.8','8.9','9.0','9.1',$
         '9.2','9.3','9.4','9.5','9.6','9.7','9.8','9.9','10.0']

  sa = n_elements(alpha)
  sY = n_elements(Y)
  stau = n_elements(tau)

  chi2_sd = dblarr(sY,stau)     ; will contain chi2 for SED
  chi2_prof = dblarr(sa)        ; will contain chi2 for 850mu, 450mu 
  chi2_sd(*,*) = !values.f_nan

  for i=0,sa-1 do begin
     dir = dusty_dir+'mod_alpha'+alpha(i)+'_y3000_tau0.1/'
     if file_test(dir+'mod_input.stb') eq 1 then begin
        p850_mod = prof_dusty(dir, scale, p850, l850, ps850, res850, w850)
        p450_mod = prof_dusty(dir, scale, p450, l450, ps450, res450, w450)
        chi2_prof(i) = chi2_profile(p850, p850_mod, p450, p450_mod)
     endif else chi2_prof(i) = !values.f_nan
  endfor

  alpham = alpha(where(chi2_prof eq min(chi2_prof, /nan)))

  print, "Alpha determined"
  print, alpham, min(chi2_prof,/nan)
  plot,float(alpha),chi2_prof,/ylog

  set_plot,'ps'
  device,filename=home+'/sources/'+name+'/'+name+'_alpha_chi2.eps',/enc
  plot,float(alpha),chi2_prof,/ylog,xtitle='!7a!17',ytitle='!7v!17!u2!n',$
       thick=3,xthick=3,ythick=3,charsize=1.2,charthick=3
  device,/close
  set_plot,'x'

  dir = dusty_dir+'mod_alpha'+alpham+'_y3000_tau0.1/'
  out = rascii(dir+'mod_input.out',skip=44,/double)
  r1 = out(3)*scale

  ymin = where(Y eq strcompress(string(fix(ceil(45./r1/100.)*100.)),/rem))
  if ymin eq [-1] then ymin = 0

;  for i=0,sY-1 do begin
  for i=ymin(0),min([sY-1,(ymin(0)+5)*2-5]) do begin
     for j=0,stau-1 do begin
        dir = dusty_dir+'mod_alpha'+alpham+'_y'+Y(i)+'_tau'+tau(j)+'/'
        if file_test(dir+'mod_input.stb') eq 1 then begin
           sed_mod_nobm = sed_dusty(dir, sed(*,nobm), intf)
           sed_mod_bm = dusty_beam(dir, sed(*,bm), beam, scale,Y(i))
           sed_mod = [sed_mod_nobm, sed_mod_bm]
           chi2_sd(i,j) = chi2_sed(sed, sed_mod)
        endif
     endfor
     print,'Y =',Y(i)
  endfor

  dum = min(chi2_sd, p, /nan)
  dum2 = array_indices(chi2_sd, p)
  Ym = Y(dum2(0))
  taum = tau(dum2(1))

  contour,alog10(chi2_sd),float(Y), float(tau),levels=findgen(40)*0.25,/follow
  oplot,[float(Ym)],[float(taum)],psym=1

  set_plot,'ps'
  device,filename=home+'/sources/'+name+'/'+name+'_Ytau_chi2.eps',/enc
  contour,alog10(chi2_sd),float(Y), float(tau),levels=findgen(40)*0.25,$
          /follow,c_thick=3,xthick=3,ythick=3,charsize=1.2,charthick=3,$
          xtitle='!18Y!17',ytitle='!7s!17!l100!n'
  device,/close
  set_plot,'x'

  print,"Y and tau determined"
  print, ym, ' ',taum,' ', dum

;--- plotting best-fit results

  dir = dusty_dir+'mod_alpha'+alpham+'_y'+Ym+'_tau'+taum+'/'
  sed_t(1,*) = sed_t(1,*)*1e-30*(3e10/sed(0,*)*1e4)
  plot_sed, home, dir, name, sed_t, intf, scale, ym

  p850_mod = prof_dusty(dir, scale, p850, l850, ps850, res850, w850)
  p450_mod = prof_dusty(dir, scale, p450, l450, ps450, res450, w450)
  sed_mod = sed_dusty(dir, sed, intf)

  plot_dusty, dir, home, name, p850, p850_mod, p450, p450_mod, sed, $
              sed_mod, alpham, Ym, taum, scale, D, intf, res850, res450

  env_dusty, dir, home, name, alpham, Ym, taum, scale, D

  print,(systime(1)-tid)/60.

  openw,unit,home+'/sources/'+name+'/input_'+name+'.dat',/get_lun
  printf,unit,name
  printf,unit,'========================='
  printf,unit,'Luminosity =',lum
  printf,unit,'Distance =',dist
  printf,unit,' '
  printf,unit,'Dusty dir =',dusty_dir
  printf,unit,' '
  printf,unit,'Profiles'
  printf,unit,'=========='
  printf,unit,'450 mu:'
  for i=0,n_elements(p450(0,*))-1 do printf,unit,p450(0,i),p450(1,i)
  printf,unit,' '
  printf,unit,'850 mu:'
  for i=0,n_elements(p850(0,*))-1 do printf,unit,p850(0,i),p850(1,i)
  printf,unit,' '
  printf,unit,'SED:'
  for i=0,n_elements(sed(0,*))-1 do printf,unit,sed(0,i),sed(1,i)
  close,unit
  free_lun,unit

end
