function err_ann_quad, image, xc, yc, rc, w, pos

; provide image, coordinates for centre, central radius 
; for the annulus and total width

; returns the average pixel value (normalized to annulus area)
; pixel values set equal to zero are not taken into account
;
; do annulus by quadrant arranged like this:
;
;    1 | 2
;   -------
;    4 | 3
;
; vector pos defines which quadrants should be included, 0 means not
; included, 1 means included.

  ann = 0.
  s = 0.

  reffi = rc-w/2.
  reffo = rc+w/2.
  tmp = ceil(reffo)

  quad = fltarr(4,4)
  quad(0,*) = [1,0,0,1]
  quad(1,*) = [0,1,0,1]
  quad(2,*) = [0,1,1,0]
  quad(3,*) = [1,0,1,0]

  for k=0,3 do begin
     if pos(k) eq 0 then begin
        image(max([0,(xc-(-1)^(quad(k,0)-1))-tmp*quad(k,0)]):$
              (xc+(-1)^(quad(k,1)-1))+tmp*quad(k,1), $
              max([0,(yc-(-1)^(quad(k,2)-1))-tmp*quad(k,2)]):$
              (yc+(-1)^(quad(k,3)-1))+tmp*quad(k,3)) = 0.
     endif
  endfor

;  for i = max([0,xc-tmp]), xc+tmp do begin
;     for j = max([0,yc-tmp]), yc+tmp do begin
  for i = xc-tmp, xc+tmp do begin
     for j = yc-tmp, yc+tmp do begin
        frac = (pixwt(xc,yc,reffo,i,j)-pixwt(xc,yc,reffi,i,j))
        ann = ann + frac*image(i,j)^2
        if image(i,j) ne 0. then s = s + frac
     endfor
  endfor

  return, sqrt(ann / s)

end
