pro image_scuba, in_name, out_name, xc, yc, rad, drad, pos, flag

  home = '/data/kristensen/wish/dusty_model/'

  if file_test(home+'scuba_data/'+in_name+'.emi.fits') eq 1 then begin

     emm = readfits(home+'scuba_data/'+in_name+'.emi.fits')

     if flag eq 450 then ps = 3 else ps = 6

     bound = 120/ps

     x = findgen(2*bound+1)*ps-120
     delta = emm(xc,yc)/10.
     levs = findgen(10)*delta

     ro = (max(rad)+drad/2.)*ps

     xb = [-1,1,1,-1]*ro/2.
     yb = [1,1,-1,-1]*ro/2

     t=strarr(3)
     t(0) = delta
     t(1) = flag
     t(2) = delta*10.
     t=strcompress(t)
     tit=out_name+', c_lev='+t(0)+', lambda='+t(1)+'mu, max='+t(2)

     set_plot,'ps'
     device,filename=home+'images/im_'+out_name+'.ps',xsize=15,ysize=15

     contour, emm(xc-bound:xc+bound,yc-bound:yc+bound), x, x, $
              levels = levs, /xst,/yst,$
              xtitle='Offset ["]', ytitle='Offset ["]', title=tit,$
              pos = [0.12,0.08,0.95,0.91]

     tvcircle, ro, 0, 0, linestyle=1, /data
     for i=0,3 do begin
        if pos(i) eq 1 then tvbox, ro, xb(i), yb(i), /data
     endfor

     device,/close
     set_plot,'x'

  endif

end
