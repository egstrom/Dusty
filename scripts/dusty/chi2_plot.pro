pro chi2_plot, chi2, pl

; assume that the chi2 matrix is the one provided by result_dusty,
; i.e. it is ordered as [[Y], [alpha], [tau], [chi2(450), chi2(850),
; chi2(SED), chi2(tot), chi2(prof)]]

  pm = !p.multi
  !p.multi = [0,2,2]

  sz = size(chi2)

;  pl = [-1,-1,12]
; which contours to plot. Use indices. '-1' indicates a variable.

  Yd = findgen(sz(1))*100+100
  alpha = [findgen(sz(2))*0.1+0.5]
  tau = findgen(sz(3))*0.1 + 0.2

  tit = ['450mu','850mu','SED','Total','Profile']
  Yt = 'Y'
  at = 'alpha'
  tt = 'tau'

  if pl(0) ne (-1) then begin
     temp = fltarr(sz(2),sz(3),sz(4))
     temp(*,*,*) = chi2(pl(0),*,*,*)
     x = alpha
     y = tau
     xt = at
     yt = tt
  endif

  if pl(1) ne (-1) then begin
     temp = fltarr(sz(1),sz(3),sz(4))
     temp(*,*,*) = chi2(*,pl(1),*,*)
     x = Yd
     y = tau
     xt = Yt
     yt = tt
  endif

  if pl(2) ne (-1) then begin
     temp = fltarr(sz(1),sz(2),sz(4))
     temp(*,*,*) = chi2(*,*,pl(2),*)
     x = Yd
     y = alpha
     xt = Yt
     yt = at
  endif

  for i=0,sz(4)-1 do begin
     contour,temp(*,*,i),x,y,levels=findgen(10)+1,$
             title=tit(i),/xst,/yst,xtitle=xt,ytitle=yt,$
             c_annotation = string(indgen(10)+1)
  endfor

  !p.multi = pm

end
