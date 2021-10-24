
FUNCTION ch_dem_gauss_fit_fn_log, x, p

;+
; NAME:
;      CH_DEM_GAUSS_FIT_FN_LOG
;
; PURPOSE:
;      A fit function in the format used by the MPFIT procedures. To
;      be used with ch_dem_gauss_fit.pro.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM; fitting.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_GAUSS_FIT_FN_LOG( X, P )
;
; INPUTS:
;      X:    Should be a simple index array that identifies the lines
;            to be computed. The line data itself is passed in a
;            common block (see below).
;      P:    A 3-element array defining the Gaussian fit function.
;              p[0]: EM_0 (amplitude)
;              p[1]: log T_0  (centroid)
;              p[2]: sigma_logT  (Gaussian width)
;
; OUTPUTS:
;      Computes line intensities (units: erg/cm2/s/sr) for a set of
;      N(X) emission lines. Note that the emission line contribution
;      functions are passed in a common block.
;
; COMMON BLOCKS:
;      CONTRIB:  Contains the contribution functions (DATA) and a log
;                temperature array (LTEMP).
;
; MODIFICATION HISTORY:
;      Ver.1, 23-Jul-2018, Peter Young.
;         Modified from ch_dem_gauss_fit_fn to use logT as the X-axis
;         rather than T.
;-

;COMMON contrib, data, ltemp, scl_aa

data=other.data
ltemp=other.ltemp
scl_aa=other.scl
ab_ref=other.ab_ref
ab_ind=other.ab_ind

dlogt=ltemp[1]-ltemp[0]

;
; Compute DEM (phi). It is defined in terms of logT rather than T.
;
phi=p[0] * scl_aa[0] * exp( -(ltemp-p[1]*scl_aa[1])^2/2./(p[2]*scl_aa[2])^2 ) / (p[2]*scl_aa[2]) / sqrt(2.*!pi)

;
; Compared to ch_dem_gauss_fit_fn, the factor ln(10)*T is missing from
; the integral.
;
s=size(data,/dim)
int=dblarr(s[0])
FOR i=0,s[0]-1 DO BEGIN
  IF n_elements(ab_ind) EQ 1 THEN BEGIN
    abfac=1.0
  ENDIF ELSE BEGIN
    abfac=p[ab_ind[i]]
  ENDELSE 
  int[i]=total(phi * abfac * ab_ref * data[*,i]  * dlogt  )
ENDFOR 

return,int

END
