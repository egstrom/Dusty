pro plot_profiles

  loadct,39

  dusty_dir = '/data1/kristensen/home/scuba/dusty/results/'

  m1 = rascii(dusty_dir+'mod_alpha1.8_y2600_tau5.8/mod_input.rtb')
  m2 = rascii(dusty_dir+'mod_alpha1.8_y2600_tau6.3/mod_input.rtb')
  m3 = rascii(dusty_dir+'mod_alpha1.8_y2600_tau6.4/mod_input.rtb')

  rin = [10.21,12.14,12.60]
  nin = [9.02d9,8.24d9,8.07d9]

  set_plot,ps
  device,filename='n1333-i4a_tprof.eps',/color,bits_per_pixel=24,/enc

  plot,[0],[0],xrange=[10,1e4],yrange=[10,300],/xst,/yst,/xlog,$
       xthick=3,ythick=3,charthick=3,charsize=1.5,$
       xtitle='!18r!17 (AU)',ytitle='!18T!17!ldust!n (K)'

  oplot,m1(0,*)*rin(0),m1(7,*),col=80,thick=5
  oplot,m2(0,*)*rin(1),m2(7,*),col=140,thick=5
  oplot,m3(0,*)*rin(2),m3(7,*),col=2540,thick=5

  device,/close
  device,filename='n1333-i4a_nprof.eps',/color,bits_per_pixel=24,/enc

  plot,[0],[0],xrange=[10,1e4],yrange=[1e4,1e10],/xst,/yst,/xlog,/ylog,$
       xthick=3,ythick=3,charthick=3,charsize=1.5,$
       xtitle='!18r!17 (AU)',ytitle='!18n!17 (cm!u-3!n)'

  oplot,m1(0,*)*rin(0),nin(0)*m1(0,*)^(-1.8),thick=5,col=80
  oplot,m2(0,*)*rin(1),nin(1)*m2(0,*)^(-1.8),thick=5,col=140
  oplot,m3(0,*)*rin(2),nin(2)*m3(0,*)^(-1.8),thick=5,col=254

  device,/close

end
