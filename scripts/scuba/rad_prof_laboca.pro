pro rad_prof_laboca, in_name, out_name, xc, yc, r0, rad, drad, pos, xn, yn, flag

  home = '/data/kristensen/wish/dusty_model/hh46/'

  emm = readfits(home+in_name+'.fits')

  emm(where(emm lt 0.))=0.

  mp = max(emm(xc-r0:xc+r0,yc-r0:yc+r0), p)

  xc = (p mod (2*r0+1)) + xc - r0
  yc = (p/(2*r0+1)) + yc - r0

  xn = xc
  yn = yc

  ps = 4.5
  fwhm = 1.38

  srad=n_elements(rad)

  tot=fltarr(srad)
  er=fltarr(srad)

  for r=0,srad-1 do tot(r) = ann_quad(emm,xc,yc,rad(r),drad,pos)

  tot = tot/tot(0)
  er(*) = 0.1*tot

  openw, unit, home+'/prof_'+out_name+'.txt', /get_lun

  printf,unit,'      b [arcsec]  I/I0         err(I/I0)'
  printf,unit,'c---------------------------------------'

  for r=0,n_elements(rad)-1 do printf,unit,[rad(r),tot(r),er(r)]
print,tot
  close, unit
  free_lun,unit

  x=findgen(150)*0.1

;  set_plot,'ps'
;  device,filename=home+'plot_'+out_name+'.ps'

  ploterror,rad*ps,tot,er,psym=1,/xlog,/ylog,xrange=[1,60],yrange=[1e-3,2],$
            /xst,/yst, title = out_name,xtitle='b [pix]',ytitle='I(b)/I(0)'

  oplot,x*ps,exp(-0.5*(x/fwhm)^2)/exp(-0.5*(rad(0)/fwhm)^2),linestyle=1

;  device,/close
;  set_plot,'x'

end
