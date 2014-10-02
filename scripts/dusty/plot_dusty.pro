pro plot_dusty, dir, home, name, p850, p850_mod, p450, p450_mod, $
                sed, sed_mod, alpham, Ym, $
                taum, scale, D, intf, res850, res450

  out = rascii(dir+'mod_input.out',skip=44,/double)
  spec = rascii(dir+'mod_input.i001',/double)
  sed_mod_t = rascii(dir+'mod_input.stb',/double)

  spec(*,0) = 0.

  r1 = out(3)*scale
  b = spec(0,*)*r1(0)

  pm_tmp = !p.multi
  !p.multi = [0,2,2]

  set_plot,'ps'
  device,filename=home+'sources/'+name+'/bf_'+name+'.eps',/enc

  ploterror,p450(0,*),p450(1,*),p450(2,*),psym=1,/xlog,/ylog,$
            xrange=[1,60],/xst,yrange=[1e-3,2],/yst,$
            xtitle='b [arcsec]', ytitle='I(b) / I(0)', $
            title='!7k!17 = 450 !7l!17m'
  oplot,p450(0,*),p450_mod(*)
  psf=0.88*exp(-.5*(b/(res450/(2.*sqrt(2.*alog(2.)))))^2) + $
      0.12*exp(-.5*(b/(40./(2.*sqrt(2.*alog(2.)))))^2)
  oplot,b,psf,linestyle=2

  ploterror,p850(0,*),p850(1,*),p850(2,*),psym=1,/xlog,/ylog,$
            xrange=[1,60],/xst,yrange=[1e-3,2],/yst,$
            xtitle='b [arcsec]', ytitle='I(b) / I(0)',$ 
            title='!7k!17 = 850 !7l!17m'
  oplot,p850(0,*),p850_mod(*)
  psf=0.88*exp(-.5*(b/(res850/(2.*sqrt(2.*alog(2.)))))^2) + $
      0.12*exp(-.5*(b/(40./(2.*sqrt(2.*alog(2.)))))^2)
  oplot,b,psf,linestyle=2

  ploterror,sed(0,*),sed(1,*),sed(2,*),/xlog,/ylog,$
            xrange=[10,1e4],psym=1,yrange=[1e-23,1e-13],/yst,$
            xtitle='!7k!17 [!7l!17m]', ytitle='!7k!17 F!i!7k!n!17 [W m!e-2!n]',$
            title='SED'
  oplot,sed_mod_t(0,*),sed_mod_t(1,*)*intf

  xyouts,[1,1,1]*1e4,[3,4,5]*1e3,['!7a!17 ='+alpham,'Y ='+Ym,'!7s!17 ='+taum],$
         /device

  device,/close
  set_plot,'x'

  !p.multi = pm_tmp

end
