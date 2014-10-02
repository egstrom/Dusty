pro plot_sed, home, dir, name, sed, intf, scale, Y

  model = rascii(dir+'mod_input.stb')

  bm = where(sed(2,*) ne 0.)

  set_plot,'ps'
  device,filename=home+'/sources/'+name+'/'+name+'_sed.eps',/enc

  plot,sed(0,*),sed(1,*),psym=1,xthick=3,ythick=3,thick=3,$
       charthick=3,charsize=1.3,xtitle='!17Wavelength (micron)',$
       ytitle='Flux (W cm!u-2!n)',/xlog,/ylog,yrange=[1e-19,1e-15],$
       xrange=[30,3e3],/xst

  sed_mod_bm = dusty_beam(dir, sed(*,bm), sed(2,bm), scale, Y)

  oplot,[sed(0,bm)],[sed_mod_bm],psym=4,thick=3

  oplot,model(0,*),model(1,*)*intf,thick=3
  oplot,model(0,*),model(1,*)*intf*2.,thick=2,col=80
  oplot,model(0,*),model(1,*)*intf/2.,thick=2,col=80

  device,/close
  set_plot,'x'

end
