function fft_expand, image, mag, filter=filter
;+
;  Name: FFT_EXPAND
;
;  Purpose: Expansion of an image with subpixelization using
;           the Fourier transform method
;           REF: Potts et al. 2003. Solar Physics, 217, 69
;  Calling Sequence:
;           Expanded Image = FFT_EXPAND (Image, Mag, filter=Filter)
;  Inputs:
;           Image             image to be expanded
;           Mag               magnification factor
;  Optional Ketword Input:
;           Filter            if set and defined, a Fourier filtering is applied.
;  History:
;           2004 May  Jongchul Chae
;-
s=size(image)
nx=s(1)
ny=s(2)

m = median(image)

if n_elements(filter) eq 0 then filter=1.
x=findgen(nx) # replicate(1., ny)
y=replicate(1., nx) # findgen(ny)
x1=x/nx
y1=y/ny

w = ( 1-((x1-0.5)/0.53)^2)*( 1- ((y1-0.5)/0.53)^2 ) ; Welch Window


f = filter*fft( w*(image-m), -1)

image1 = fltarr(nx*mag, ny*mag, /nozero)
kx = 2*!pi*[findgen(nx/2+1),reverse(-1-findgen(nx-nx/2-1))]/nx
ky = 2*!pi*[findgen(ny/2+1),reverse(-1-findgen(ny-ny/2-1))]/ny
kx=kx#replicate(1, ny)
ky=replicate(1,nx)#ky


for i = 0, mag-1 do for j=0, mag-1 do begin

dx  = float(i)/mag
dy  = float(j)/mag

fs = float(fft(f * exp(complex(0., (kx*dx+ky*dy))), 1))
x1=(x+dx)/nx
x1= x1*(x1 ge 0 and x1 le 1) + (x1-1)*(x1 gt 1.)+(x1+1)*(x1 lt 0)
y1=(y+dx)/ny
y1= y1*(y1 ge 0 and y1 le 1) + (y1-1)*(y1 gt 1)+(y1+1)*(y1 lt 0)
ws = ( 1-((x1-0.5)/0.53)^2)*( 1- ((y1-0.5)/0.53)^2)



image1(x*mag+i, y*mag+j) = fs/ws + m

endfor

return, image1

end