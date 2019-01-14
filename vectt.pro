PRO vectt,UU,VV,dat=dat,ISARROW,SX,SY,X,Y, Missing = Missing, Length = length, $
      Dots = dots,Thick=thick, Rlen=rlen, T3d=t3d, $
       _EXTRA = extra
;
;+
; NAME:
;	VELOVECT
;
; PURPOSE:
;	Produce a two-dimensional velocity field plot.
;
;	A directed arrow is drawn at each point showing the direction and
;	magnitude of the field.
;
; CATEGORY:
;	Plotting, two-dimensional.
;
; CALLING SEQUENCE:
;	VELOVECT, U, V [, X, Y]
;
; INPUTS:
;	U:	The X component of the two-dimensional field.
;		U must be a two-dimensional array.
;
;	V:	The Y component of the two dimensional field.  Y must have
;		the same dimensions as X.  The vector at point (i,j) has a
;		magnitude of:
;
;			(U(i,j)^2 + V(i,j)^2)^0.5
;
;		and a direction of:
;
;			ATAN2(V(i,j),U(i,j)).
;
; OPTIONAL INPUT PARAMETERS:
; 	X:	Optional abcissae values.  X must be a vector with a length
;		equal to the first dimension of U and V.
;
;	Y:	Optional ordinate values.  Y must be a vector with a length
;		equal to the first dimension of U and V.
;
; KEYWORD INPUT PARAMETERS:
;      MISSING:	Missing data value.  Vectors with a LENGTH greater
;		than MISSING are ignored.
;
;	LENGTH:	Length factor.  The default of 1.0 makes the longest (U,V)
;		vector the length of a cell.
;
;	DOTS:	Set this keyword to 1 to place a dot at each missing point.
;		Set this keyword to 0 or omit it to draw nothing for missing
;		points.  Has effect only if MISSING is specified.
;
;	COLOR:	The color index used for the plot.
;
;       RLEN:   The length of arrow head. If rlen=0.0, then no arrow.
;
;	Note:   All other keywords are passed directly to the PLOT procedure
;		and may be used to set option such as TITLE, POSITION,
;		NOERASE, etc.
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Plotting on the selected device is performed.  System
;	variables concerning plotting are changed.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Straightforward.  Unrecognized keywords are passed to the PLOT
;	procedure.
;
; MODIFICATION HISTORY:
;	DMS, RSI, Oct., 1983.
;	For Sun, DMS, RSI, April, 1989.
;	Added TITLE, Oct, 1990.
;	Added POSITION, NOERASE, COLOR, Feb 91, RES.
;	August, 1993.  Vince Patrick, Adv. Visualization Lab, U. of Maryland,
;		fixed errors in math.
;	August, 1993. DMS, Added _EXTRA keyword inheritance.
;-
;
      on_error,2                      ;Return to caller if an error occurs
	  u=uu
	  v=vv
        s = size(uu)
        t = size(vv)
        if s(0) ne 2 then begin
baduv:   message, 'U and V parameters must be 2D and same size.'
                endif
        if total(abs(s(0:2)-t(0:2))) ne 0 then goto,baduv
;
        PRINT,'PARAMETER NUMBERS:',N_PARAMS()
        if n_params()  lt 3  then isarrow=1.
        if n_params()  lt 4  then  SX=S(1)
        if n_params()  lt 5  then  SY=S(2)
        if n_params(0) lt 6 then x = findgen(SX) else $
                if n_elements(x) ne s(1) then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 7 then y = findgen(SY) else $
                if n_elements(y) ne s(2) then goto,badxy

;
        if n_elements(missing) le 0 then missing = 1.0e30
        if n_elements(length) le 0 then length = 1.0
        if n_elements(rlen) le 0 then   rlen=isarrow*0.3
	if n_elements(t3d) le 0 then t3d=0
	if n_elements(thick) le 0 then thick=1

        IF (SX NE S(1))  OR  (SY NE S(2)) THEN BEGIN
            U=CONGRID(UU,SX,SY,_extra=extra)
            V=CONGRID(VV,SX,SY,_extra=extra)
	    dat=congrid(dat,sx,sy)
        ENDIF


        mag = sqrt(u^2+v^2)             ;magnitude.
                ;Subscripts of good elements
        nbad = 0                        ;# of missing points
        if n_elements(missing) gt 0 then begin
                good = where(mag lt missing)
                if keyword_set(dots) then bad = where(mag ge missing, nbad)
        endif else begin
                good = lindgen(n_elements(mag))
        endelse

        ugood = u(good)
        vgood = v(good)
        x0 = min(x)                     ;get scaling
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)
	x_step=(x1-x0)/sx
	y_step=(y1-y0)/sy
	maxmag=max([max(ugood/x_step),max(vgood/y_step)])
	sina = length * (ugood/maxmag)
	cosa = length * (vgood/maxmag)
;
        if n_elements(title) le 0 then title = ''
        ;--------------  plot to get axes  ---------------
        if n_elements(color) eq 0 then color = !p.color
        x_b0=x0-x_step
	x_b1=x1+x_step
	y_b0=y0-y_step
	y_b1=y1+y_step
        if n_elements(position) eq 0 then begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            color=color, t3d=t3d, _EXTRA = extra
        endif else begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            color=color, t3d=t3d, _EXTRA = extra
        endelse
;
;        if keyword_set(arrow) then r=0.0 else r=.3                                       ;        r = .3                           ;len of arrow head
         r = rlen
        angle = 22.5 * !dtor            ;Angle of arrowhead
        st = r * sin(angle)             ;sin 22.5 degs * length of head
        ct = r * cos(angle)
;
        for i=0L,n_elements(good)-1 do begin     ;Each point
                x0 = x(good(i) mod sx)        ;get coords of start & end
                dx = sina(i)
                x1 = x0 + dx
                y0 = y(good(i) / sx)
                dy = cosa(i)
                y1 = y0 + dy
		xd=x_step
		yd=y_step
;              ---------------------------
                if dat(x0,y0) gt 0.0 then color0=0
                if dat(x0,y0) lt 0.0 then color0=255
                if dat(x0,y0) eq 0.0 then color0=127
                if dx eq 0 and dy eq 0 then goto,label
;              ----------------------------
                plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
			x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
                      [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
			y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
                      color=color0,t3d=t3d, thick=thick
label:
                endfor
;        if nbad gt 0 then $             ;Dots for missing?
;                oplot, x(bad mod sx), y(bad / sx), psym=3, $
;                    color=color, t3d=t3d, thick=thick
end
