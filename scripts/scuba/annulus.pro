function annulus, image, xc, yc, rc, w

; provide image, coordinates for centre, central radius 
; for the annulus and total width

; returns the average pixel value (normalized to annulus area)
; pixel values set equal to zero are not taken into account

ann = 0.
s = 0.

reffi = rc-w/2.
reffo = rc+w/2.
tmp = ceil(reffo)

for i = xc-tmp, xc+tmp do begin
   for j = yc-tmp, yc+tmp do begin
      frac = (pixwt(xc,yc,reffo,i,j)-pixwt(xc,yc,reffi,i,j))
      ann = ann + frac*image(i,j)
      if image(i,j) ne 0. then s = s + frac
   endfor
endfor

return, ann / s

end
