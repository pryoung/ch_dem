
FUNCTION ch_dem_linear_fit_fn, x, p, _extra=e


;+
; NAME:
;      CH_DEM_LINEAR_FIT_FN
;
; PURPOSE:
;      A fit function in the format used by the MPFIT procedures. To
;      be used with ch_dem_linear_fit.pro.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM; fitting.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_LINEAR_FIT_FN( X, P )
;
; INPUTS:
;      X:    Should be a simple index array that identifies the lines
;            to be computed. 
;      P:    A 1D array containing the fit parameters.
;      _EXTRA: This is expected to be a structure containing auxiliary
;              information needed to compute the intensities. It is
;              defined in ch_dem_linear_fit. 
;
; OUTPUTS:
;      Computes line intensities (units: erg/cm2/s/sr) for a set of
;      N(X) emission lines. 
;
; MODIFICATION HISTORY:
;      Ver.1, 18-Jul-2019, Peter Young
;         I used ch_dem_gauss_fit_fn as a starting point.
;-




phi=ch_dem_linear_compute_dem(e.ltemp,e.ltmp_nodes,p,scl=e.scl,em=e.em)


;
; I compute the intensities from the DEM using
; ch_dem_compute_ints. This requires a structure containing various
; auxiliary information.
;
data={ phi: phi, $
       line_data: e.line_data, $
       abstr: e.abstr, $
       ab_ref: e.ab_ref, $
       ltemp: e.ltemp, $
       p: p }
;      
int=ch_dem_compute_ints(data)


return,int

END

