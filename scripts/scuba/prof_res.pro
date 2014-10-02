pro prof_res

  home = '/data/kristensen/wish/dusty_model/'

  data = rascii(home+'scuba_position_supp.txt', header = in_name)
  temp = rascii(home+'name_supp.txt', header = out_name)

  stot = n_elements(in_name)

  rad450 = [2,4,6,8,10,12,14,16,18]
  rad850 = [1.5,3,4.5,6,7.5,9]
  drad450 = 2.
  drad850 = 1.5

  for i=0,stot-1 do begin
     in=in_name(i)+'_450um'
     flag=450
     rad_prof_scuba, in, out_name(i)+'_450mu', 2*data(0,i), 2*data(1,i), $
                     1, rad450, drad450, data(2:5,i), xn, yn, flag
     image_scuba, in, out_name(i)+'_450mu', xn, yn, $
                  rad450, drad450, data(2:5,i), flag

     in=in_name(i)+'_850um'
     flag=850
     rad_prof_scuba, in, out_name(i)+'_850mu', data(0,i), data(1,i), $
                     1, rad850, drad850, data(2:5,i), xn, yn, flag
     image_scuba, in, out_name(i)+'_850mu', xn, yn, $
                  rad850, drad850, data(2:5,i), flag
  endfor

end
