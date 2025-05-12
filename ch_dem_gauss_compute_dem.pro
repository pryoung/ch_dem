
FUNCTION ch_dem_gauss_compute_dem, ltemp, p, scl=scl, temp_log=temp_log


;+
; NAME:
;     CH_DEM_GAUSS_COMPUTE_DEM
;
; PURPOSE:
;     Computes the DEM values for a set of input temperatures. This is a
;     Gaussian DEM, defined by three parameters.
;
; CATEGORY:
;     CHIANTI; differential emission measure; DEM; fitting.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_GAUSS_COMPUTE_DEM( Ltemp, P )
;
; INPUTS:
;     Ltemp:  An array of log temperatures.
;     P:      A 3-element array containing the fit parameters for the
;             Gaussian DEM.
;
; OPTIONAL INPUTS:
;     Scl:    A 3-element array containing scaling parameters such that
;             SCL*P are the true DEM parameters.
;	
; KEYWORD PARAMETERS:
;     TEMP_LOG: If set, then the DEM is defined in logT space rather than T.
;
; OUTPUTS:
;     The DEM function defined for the input temperatures LTEMP.
;
; MODIFICATION HISTORY:
;     Ver.1, 12-May-2025, Peter Young
;-


IF n_elements(scl) EQ 0 THEN scl=fltarr(3)+1.

t=10.^ltemp

;
; Compute DEM (phi). This expression comes from Eq. 3 of Warren &
; Brooks (2009, ApJ, 700, 762). The temp_log function comes from
; Phi(T)*dT/d(logT).
;
IF keyword_set(temp_log) THEN BEGIN
  phi=p[0] * scl[0] * exp( -(ltemp-p[1]*scl[1])^2/2./(p[2]*scl[2])^2 ) / (p[2]*scl[2]) / sqrt(2.*!pi) $
      / (t* alog(10.) )
ENDIF ELSE BEGIN
  phi=p[0] * scl[0] * exp( -(t-p[1]*scl[1])^2/2./(p[2]*scl[2])^2 ) / (p[2]*scl[2]) / sqrt(2.*!pi)
ENDELSE

return,phi

END
