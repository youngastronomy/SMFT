pro crossfct, img1, img2, x, y, mx, crfct
;+
; NAME:
;       CROSSFCT
;
; PURPOSE:
;       This procedure computes shifts between two images by means of
;       the crosscorrelation function.
;
; CATEGORY:
;       Image processing.
;
; CALLING SEQUENCE:
;       CROSSFCT, img1, img2, x, y [, mx, crfct]
;
; INPUTS:
;       IMG1:    Two-dimensional image.
;
;       IMG2:    Two-dimensional image. The images must have the same
;                dimensions and the dimensions have to be even numbers.
;                Since the procedure makes use of Fourier transforms,
;                it is recommended to use dimensions which are multiples
;                of small prime numbers.
;
; OUTPUTS:
;       X:       Shifts in x-direction with respect to IMG1
;       Y:       Shifts in y-direction with respect to IMG1
;       MX:      Maximimum of the crosscorrelation function
;       CRFCT:   Crosscorrelation function of IMG1 and IMG2
;
;
; MODIFICATION HISTORY:
;       04-Sep-1997    C. Denker     BBSO/NJIT
;-

   ;;; Get the image size.
dim1 = SIZE ( img1 )
dim2 = SIZE ( img2 )
nx = dim1 ( 1 )
ny = dim1 ( 2 )

   ;;; Compute the crosscorrelation function.
crfct = FFT ( img1, -1 )
crfct = CONJ ( crfct )
crfct = crfct * FFT ( img2, -1 )
crfct = FFT (crfct, 1 )
crfct = FLOAT ( crfct )

   ;;; Compute the shifts with respect to IMG1.
mx = MAX ( crfct, index )
x = index MOD nx
IF x GT ( nx / 2 ) THEN x = x - nx
y = index / nx
IF y GT ( ny / 2 ) THEN y = y - ny

   ;;; Shift crosscorelation function from butterfly to centered format.
crfct = SHIFT ( crfct, nx / 2, ny / 2 )

RETURN

END

;==================================================================================
