pro rad_prof_scuba, in_name, out_name, xc, yc, r0, rad, drad, pos, xn, yn, flag

  home = '/data/kristensen/wish/dusty_model/'

  if file_test(home+'scuba_data/'+in_name+'.emi.fits') eq 1 then begin

     emm = readfits(home+'scuba_data/'+in_name+'.emi.fits')
     err = readfits(home+'scuba_data/'+in_name+'.err.fits')

     emm(where(emm lt 0.))=0.
     err(where(emm eq 0.))=0.

     mp = max(emm(xc-r0:xc+r0,yc-r0:yc+r0), p)

     xc = (p mod (2*r0+1)) + xc - r0
     yc = (p/(2*r0+1)) + yc - r0

     xn = xc
     yn = yc

     if flag eq 450 then ps = 3 else ps = 6
     if flag eq 450 then fwhm = 2.12 else fwhm = 1.38

     srad=n_elements(rad)

     tot=fltarr(srad)
     er=fltarr(srad)

     for r=0,srad-1 do begin
        tot(r) = ann_quad(emm,xc,yc,rad(r),drad,pos)
        er(r) = err_ann_quad(err,xc,yc,rad(r),drad,pos)
     endfor

     er = er/tot(0)
     tot = tot/tot(0)

     openw, unit, home+'profiles/prof_'+out_name+'.txt', /get_lun

     printf,unit,'      b [arcsec]  I/I0         err(I/I0)'
     printf,unit,'c---------------------------------------'

     for r=0,n_elements(rad)-1 do printf,unit,[rad(r),tot(r),er(r)]

     close, unit
     free_lun,unit

     x=findgen(150)*0.1

     set_plot,'ps'
     device,filename=home+'profiles/plot_'+out_name+'.ps'

     ploterror,rad*ps,tot,er,psym=1,/xlog,/ylog,xrange=[1,60],yrange=[1e-3,2],$
               /xst,/yst, title = out_name,xtitle='b [pix]',ytitle='I(b)/I(0)'

     oplot,x*ps,exp(-0.5*(x/fwhm)^2)/exp(-0.5*(rad(0)/fwhm)^2),linestyle=1

     device,/close
     set_plot,'x'

  endif

end
