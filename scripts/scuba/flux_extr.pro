function flux_extr, in_name, xc, yc, r0, rad, pos

  home = '/data/kristensen/wish/dusty_model/'

  if file_test(home+'scuba_data/'+in_name+'.emi.fits') eq 1 then begin

     emm = readfits(home+'scuba_data/'+in_name+'.emi.fits')

     emm(where(emm lt 0.))=0.

     mp = max(emm(xc-r0:xc+r0,yc-r0:yc+r0), p)

     xc = (p mod (2*r0+1)) + xc - r0
     yc = (p/(2*r0+1)) + yc - r0

     xn = xc
     yn = yc

     tot = ann_quad(emm,xc,yc,rad/2.,rad,pos)

     return, tot

  endif else begin
     return, !values.f_nan
  endelse


end
