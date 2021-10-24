
FUNCTION ch_dem_gauss_fit_fn, x, p, _extra=e

;+
; NAME:
;      CH_DEM_GAUSS_FIT_FN
;
; PURPOSE:
;      A fit function in the format used by the MPFIT procedures. To
;      be used with ch_dem_gauss_fit.pro.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM; fitting.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_GAUSS_FIT_FN( X, P )
;
; INPUTS:
;      X:    Should be a simple index array that identifies the lines
;            to be computed. 
;      P:     A 1D array containing the fit parameters. The first
;             three parameters define the Gaussian (amplitude,
;             centroid, width). Remaining parameters are for element
;             abundances. 
;      _EXTRA: This is expected to be a structure containing auxiliary
;              information needed to compute the intensities. It is
;              defined in ch_dem_gauss_fit. 
;
; KEYWORD PARAMETERS:
;      TEMP_LOG: If set, then the Gaussian is defined in log-T space
;                rather than T.
;
; OUTPUTS:
;      Computes line intensities (units: erg/cm2/s/sr) for a set of
;      N(X) emission lines. 
;
; MODIFICATION HISTORY:
;      Ver.1, 18-Jul-2018, Peter Young.
;      Ver.2, 19-Jul-2018, Peter Young.
;         Added scaling parameter (scl_aa) to the common block.
;      Ver.3, 8-Feb-2019, Peter Young.
;         I've removed the common block and now pass the extra
;         data through a structure called OTHER; the function now
;         allows for the element abundances to be varied as part of
;         the minimization process.
;      Ver.4, 18-Jul-2019, Peter Young
;         Replaced OTHER keyword with _EXTRA; updated header.
;      Ver.5, 29-Jul-2019, Peter Young
;         Moved "extra_fac" outside of the nl intensity loop and
;         re-defined phi instead; removed /log keyword (now passed in
;         the _extra structure); now ch_dem_compute_ints to compute
;         intensities. 
;-



phi=ch_dem_gauss_compute_dem(e.ltemp,p,scl=e.scl,temp_log=e.temp_log)


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

