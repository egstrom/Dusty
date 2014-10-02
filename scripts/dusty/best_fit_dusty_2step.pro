pro best_fit_dusty_2step, home, name, chi2, alpha, Y, tau, alpham, Ym, taum

; determine location of chi2 minima and give confidence intervals for
; best-fit parameters

  schi2 = size(chi2)
  sY = schi2(1)
  sa = schi2(2)
  stau = schi2(3)

  Yno = fltarr(sY)
  ano = fltarr(sa)
  tno = fltarr(stau)

  Yno(*) = Y(*)
  ano(*) = alpha(*)
  tno(*) = tau(*)

  title = ['450mu', '850mu', '  SED', 'Total', 'Profile']

  openw,unit,home+'results/bf_sum_'+name+'.txt',/get_lun

  t = min(chi2(*,*,*,4), p, /nan)

  dum = array_indices(chi2(*,*,*,4), p)

  alpham = alpha(dum(1))        ; determine slope first, then fix it 
;   loc_alpha=where(alpha eq alpham)

  t = min(chi2(*,dum(1),*,2), p, /nan)
  dum2 = array_indices(chi2(*,dum(1),*,2), p)
print,dum2(0),dum(1),dum2(2)
  Ym = Y(dum2(0))
  taum = tau(dum2(2))

  conf = where(chi2(*,*,*,4) le t+1.)

  if conf eq [-1] then begin
     print,'Model caused problems'
     printf,unit,'-----------------------------------------------'
     printf, unit, '!!! Singularity detected !!!'
     printf, unit, 'Best fit parameters are'
     printf,unit,['  chi2(min) = ', string(format='(f6.2)',t)]
     for j=0, schi2(4)-1 do begin
        if j ne 3 then printf,unit, ['chi2('+title(j)+') = ', $
                                     string(format='(f6.2)', $
                                            chi2(dum2(0), dum(1), $
                                                 dum2(2),j))]
     endfor

     printf,unit, ['        Y = ', string(format='(f6.0)',Ym)]
     printf,unit, ['    alpha = ',string(format='(f6.1)',alpham)]
     printf,unit, ['      tau = ', string(format='(f6.1)',taum)]

  endif else begin

     printf,unit,'-----------------------------------------------'
     printf,unit,title(3)
     printf,unit,['  chi2(min) = ', string(format='(f6.2)',t)]
     for j=0, 2 do begin
        printf,unit, ['chi2('+title(j)+') = ', $
                      string(format='(f6.2)', $
                             chi2(dum2(0), dum(1), dum2(2),j))]
     endfor

     printf,unit, ['        Y = ', string(format='(f6.0)',Ym), $
                   '  [',string(format='(f6.0)', min(Yno(conf mod sY), /nan)), '; ',$
                   string(format='(f6.0)',max(Yno(conf mod sY), /nan)), ']']
     printf,unit, ['    alpha = ',string(format='(f6.1)',alpham) $
                   , '  [', string(format='(f6.1)', min(ano((conf mod (sY*sa))/sY), /nan)) $
                   , '; ', string(format='(f6.1)', max(ano((conf mod (sY*sa))/sY), /nan)), ']']
     printf,unit, ['      tau = ', string(format='(f6.1)',taum), $
                   '  [', string(format='(f6.1)', min(tno(conf/(sY*sa)), /nan)), $
                   '; ', string(format='(f6.1)', max(tno(conf/(sY*sa)), /nan)), ']']
  endelse
     
  close,unit
  free_lun,unit

end
