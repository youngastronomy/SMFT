pro free_field,bz,ha,nz,pbx,pby,pbz 
ss=size(bz)
sizex=ss(1)
sizey=ss(2)
i=complex(0,1)
sub=indgen(sizex,sizey)

n=sub/sizex
m=sub-n*sizex
bfft=fft(bz,1)
k1=2.0*!pi/float(sizex)
k2=2.0*!pi/float(sizey)

kk=(sin(k1*m))^2+(sin(k2*n))^2-ha*ha
 
ssign=((kk gt 0)-0.5)*2
kk=ssign*kk
ssi=(ssign eq 1)+i*(ssign eq -1)
k=i*ssi*sqrt(kk)

zk=exp(i*k*nz)
fm=(sin(k1*m))^2+(sin(k2*n))^2+0.000001

bxfft=-zk*bfft*(-k*sin(k1*m)+i*ha*sin(k2*n))/fm
byfft=-zk*bfft*(-k*sin(k2*n)-i*ha*sin(k1*m))/fm
bzfft=zk*bfft
pbx=float(fft(bxfft,-1))
pby=float(fft(byfft,-1))
pbz=float(fft(bzfft,-1))
end
