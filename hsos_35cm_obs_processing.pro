pro hsos_35cm_obs_processing


device,decomposed=0
loadct,0


fileL=dialog_pickfile(Filter='L5*.fit')
print,fileL
fileQ=dialog_pickfile(Filter='Q5*.fit')
print,fileQ
fileU=dialog_pickfile(Filter='U5*.fit')
print,fileU

dataV=float(readfits(fileL,h))
dataQ=float(readfits(fileQ))
dataU=float(readfits(fileU))

time=sxpar(h,'STARTOBS')


U=10000.*(dataU[*,*,0]-dataU[*,*,1])/(dataU[*,*,0]+dataU[*,*,1])
Q=10000.*(dataQ[*,*,0]-dataQ[*,*,1])/(dataQ[*,*,0]+dataQ[*,*,1])
V=10000.*(dataV[*,*,0]-dataV[*,*,1])/(dataV[*,*,0]+dataV[*,*,1])



;-----------------------align data-------------------

sw=size(V)
nx=sw(1)
ny=sw(2)



window,0,xs=fix(nx*0.5),ys=fix(ny*0.5),xpos=0,ypos=0,title='L'
tvscl,congrid(dataV[*,*,0],fix(nx*0.5),fix(ny*0.5))

x0=115
y0=70
NNx=256
NNy=256
;choosing region (usually the large umbra) to align the Q,U,V each other
box_cursor,x0,y0,NNx,NNy,/message,/INIT
x0=2*x0
y0=2*Y0
NNx=2*NNx
NNy=2*NNY
crossfct,dataV[x0:x0+NNx,Y0:Y0+NNy,0]*0.5,dataU[x0:x0+NNx,Y0:Y0+NNy,0],shqx,shqy
crossfct,dataV[x0:x0+NNx,Y0:Y0+NNy,0]*0.5,dataU[x0:x0+NNx,Y0:Y0+NNy,0],shux,shuy

print,shqx,shqy
print,shux,shuy

Q=shift(Q,-shqx,-shqy)
U=shift(U,-shux,-shuy)


;--------------------cut the edge with no bad signal------------------------------
  region=[0,124,991,915]

  X0=region[0]
  Y0=region[1]
  X1=region[2]
  Y1=region[3]

Q=Q[X0:X1,Y0:Y1]
U=U[X0:X1,Y0:Y1]
V=V[X0:X1,Y0:Y1]


sw=size(V)
nx=sw(1)
ny=sw(2)



print,median(Q),median(U),median(V)

;-----------------------------try to remove the flat signal-------------------------
Q=Q-median(Q)
U=U-median(U)
V=V-median(V)





;--------------------------------show Q,U signal--------------------
window,0,xs=0.5*nx,ys=0.5*ny,title='Q',xpos=0,ypos=180
tvscl,rebin(smooth(Q,2,/edge_),0.5*nx,0.5*ny)
window,1,xs=0.5*nx,ys=0.5*ny,title='U',xpos=0.5*nx,ypos=180
tvscl,rebin(smooth(U,2,/edge_),0.5*nx,0.5*ny)
window,2,xs=0.5*nx,ys=0.5*ny,title='V',xpos=1.*nx,ypos=180
tvscl,rebin(smooth(V,2,/edge_),0.5*nx,0.5*ny)


;----------------------------smooth to improve the S/N ratio
Q=rebin(smooth(Q,2,/edge_),0.5*nx,0.5*ny)
U=rebin(smooth(U,2,/edge_),0.5*nx,0.5*ny)
V=rebin(smooth(V,2,/edge_),0.5*nx,0.5*ny)



;-----------------------------To the correct X-Y coodinates from CCD data coordinate
Q=rotate(Q,7)
U=rotate(U,7)
V=rotate(V,7)




;--------------------------------using the HSOS regular calibration

Bz0=V
l5calib,Bz0

;Q,U calibration
q5=Q
u5=U
q5=float(q5)
u5=float(u5)
bt=sqrt(sqrt(q5^2+u5^2))*93.7/sqrt(2.)
bang=!pi/2.0-0.5*atan(u5,q5)
bx0=-bt*cos(bang)
by0=bt*sin(bang)


;--------------------------Using Potential approach to remove the ambiguity----------------------
ldat=Bz0
free_field,smooth(ldat,1,/edge_),0,0,pbx,pby,pbz
test=(((pbx*bx0+pby*by0) gt 0.0)-0.5)*2
bx0=bx0*test
by0=by0*test
btt=bx0*bx0+by0*by0

;drawn the picture
sw=size(bz0)
nnx=sw(1)
nny=sw(2)


;---------------------------show the vector magnetogram---------------------------------
ratio=4
window,3,xsize=ratio*nnx,ysize=ratio*nny,XPos=75,YPos=150,title='magetic filed 35cm'


TVSCL, congrid(bz0,ratio*nnx,ratio*nny)>(-200)<200


vectt,bx0,by0,dat=fft_expand(bz0,ratio),1,60,60,position=[0,0,1.0,1.0],length=1.0,rlen=0.3,dots=0,thick=1,/noerase


jump:

help,bz0,bx0,by0

end
