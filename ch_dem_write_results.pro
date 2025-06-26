
PRO ch_dem_write_results, ld_fit, abstr, interr_scale=interr_scale


;+
; NAME:
;     CH_DEM_WRITE_RESULTS
;
; PURPOSE:
;     Prints the results of the ch_dem optimization procedures to the
;     IDL input window. Usually called from within the CH_DEM fitting
;     routine (e.g., ch_dem_linear_fit). 
;
; CATEGORY:
;     CHIANTI; DEM; output.
;
; CALLING SEQUENCE:
;     CH_DEM_WRITE_RESULTS, Ld_Fit, Abstr
;
; INPUTS:
;     Ld_Fit:  A structure in the format returned by
;              ch_dem_read_line_ids but containing only those spectral
;              features that were used for the fit.
;     Abstr:   A structure in the format returned by
;              ch_dem_process_abund containing the abundance results
;              from the DEM method.
;
; OPTIONAL INPUTS:
;     Interr_scale:  This is passed through from the DEM routines and
;                     modifies the intensity error value printed to the
;                     screen.
;
; OUTPUTS:
;     Prints text to the IDL input window giving the results from the
;     CH_DEM procedure. 
;
; MODIFICATION HISTORY:
;     Ver.1, 16-Jun-2021, Peter Young
;     Ver.2, 22-May-2025, Peter Young
;       Added int_err_scale= optional input.
;     Ver.3, 25-Jun-2025, Peter Young
;       Changed int_err_scale to interr_scale to be consistent with
;       other routines.
;     Ver.4, 26-Jun-2025, Peter Young
;       Modified how the abundances are printed to take account of
;       abund_lower and abund_upper tags (if available).
;-


;
; Print out the comparison between observed and fitted
; intensities.
;
nfit=n_elements(ld_fit)
print,'    Ion                 Wvl   Obs_Int   Obs_Err   Model_Int'
FOR i=0,nfit-1 DO BEGIN
  IF n_elements(interr_scale) NE 0 THEN BEGIN
    err=sqrt(ld_fit[i].err^2 + (ld_fit[i].int*interr_scale)^2)
  ENDIF ELSE BEGIN
    err=ld_fit[i].err
  ENDELSE 
  print,format='(a7,a20,2f10.2,f12.2)', ld_fit[i].ion,ld_fit[i].label, $
        ld_fit[i].int,err,ld_fit[i].model_int
ENDFOR

;
; Now print abundance results
;
nab=n_elements(abstr)
print,''
print,' Element     Type      Log Ab    Abundance (x10^6)'
FOR i=0,nab-1 DO BEGIN
  z2element,abstr[i].elt_num,name,/symbol
  CASE abstr[i].type OF
    0: type='algebraic'
    1: type='reference'
    2: type='optimized'
    3: type='fixed'
  ENDCASE
  IF abstr[i].type EQ 2 OR abstr[i].type EQ 0 THEN BEGIN
    IF abstr[i].error NE -1. THEN BEGIN 
      errstr=' +/- '+trim(string(format='(f6.2)',abstr[i].error*1e6))
      format='(6x,a2,a12,f8.2,f10.2,a11)'
      errstr=strpad(errstr,11,/after,fill=' ')
    ENDIF ELSE BEGIN
      errstr=' (+ '+trim(string(format='(f6.2)',abstr[i].abund_upper*1e6-abstr[i].abund*1e6))+')'+ $
             ' (- '+trim(string(format='(f6.2)',abstr[i].abund*1e6-abstr[i].abund_lower*1e6))+')'
      errstr=strpad(errstr,22,/after,fill=' ')
      format='(6x,a2,a12,f8.2,f10.2,a22)'
    ENDELSE 
  ENDIF ELSE BEGIN
    errstr='          '
  ENDELSE 
  print,format=format,name,type,alog10(abstr[i].abund)+12., $
        abstr[i].abund*1e6,errstr
ENDFOR


END
