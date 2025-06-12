
FUNCTION ch_dem_get_ltemp, line_data, dlogt=dlogt, ioneq_file=ioneq_file, $
                           log_press=log_press,log_dens=log_dens, $
                           brooks=brooks

;+
; NAME:
;     CH_DEM_GET_LTEMP
;
; PURPOSE:
;     Computes the log temperature array that is used for the DEM
;     calculation. The array is calculated based on the contribution
;     functions of the ions.
;
; CATEGORY:
;     CHIANTI; DEM.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_GET_LTEMP( Line_Data )
;
; INPUTS:
;     Line_Data:  An IDL structure in the format returned by
;                 ch_dem_read_line_ids that contains the list of emission
;                 lines and their IDs.
;
; OPTIONAL INPUTS:
;     Dlogt:  The spacing for the log temperature array. The default is
;             0.1 dex.
;     Ioneq_File:  The name of a CHIANTI ioneq file. If not specified,
;                  then it is calculated.
;     Log_dens:  Specifies the logarithm of the electron number density
;              (units: cm^-3) to be used for the calculation.Either
;              log_dens or log_press should be specified (not both).
;     Log_press:  Specifies the logarithm of the electron pressure
;              (units: K cm^-3) to be used for the calculation. Either
;              log_dens or log_press should be specified (not both).
;	
; KEYWORD PARAMETERS:
;     BROOKS:  If set, then the method of David Brooks is used for computing
;              the temperature range.
;
; OUTPUTS:
;     A 1D array containing the log temperatures. The temperatures are
;     evenly spaced in 0.1 dex intervals. The default is compute where
;     the contribution function is a factor 0.005 greater than the max
;     of the contribution function, for each line. The temperature range is
;     then min and max of all of these ranges. If the /brooks keyword is set,
;     then the David Brooks method is used whereby the temperature range of
;     each ion is set to be +/- 0.3 from log(Tmax). Ltemp is then computed
;     from the min and max of all of these ranges.
;
; EXAMPLE:
;     IDL> ltemp=ch_dem_get_ltemp(line_data,log_dens=9.0)
;     IDL> ltemp=ch_dem_get_ltemp(line_data,log_dens=9.0,/brooks)
;     IDL> ltemp=ch_dem_get_ltemp(line_data,log_press=15.0,ioneq_file=!ioneq_file)
;
; MODIFICATION HISTORY:
;     Ver.1, 11-Jun-2025, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ltemp=ch_dem_get_ltemp( line_data [,log_dens=, log_press=, /brooks'
  print,'                          ioneq_file=, dlogt= ] )'
  return,-1
ENDIF 

IF n_elements(dlogt) EQ 0 THEN dlogt=0.10

; Check the log_dens and log_press inputs. We need one or the other, but not
; both (or neither).
;
IF n_elements(log_dens) NE 0 AND n_elements(log_press) NE 0 THEN BEGIN
  message,/info,/cont,'Please specify either LOG_DENS or LOG_PRESS, but not both. Returning...'
  return,-1
ENDIF 
;
IF n_elements(log_dens) EQ 0 AND n_elements(log_press) EQ 0 THEN BEGIN
  message,/info,/cont,'Please specify either LOG_DENS or LOG_PRESS. Returning...'
  return,-1
ENDIF


ld=line_data

;
; Get contribution functions.
;
ltemp=findgen(81)/20.+4.0
ld=ch_dem_add_contrib(ld, ltemp, avalfile=avalfile, $
                      log_press=log_press, log_dens=log_dens, $
                      ioneq_file=ioneq_file, dir_lookup=dir_lookup, $
                      truncate=0)
n=n_elements(ld)

ltemp_min=8.0
ltemp_max=4.0

trunc_factor=0.005

IF keyword_set(brooks) THEN BEGIN
  FOR i=0,n-1 DO BEGIN
    contrib=ld[i].contrib
    getmax=max(contrib,imax)
    logtmax=ltemp[imax]
    IF logtmax-0.3 LT ltemp_min THEN ltemp_min=logtmax-0.3
    IF logtmax+0.3 GT ltemp_max THEN ltemp_max=logtmax+0.3
  ENDFOR
ENDIF ELSE BEGIN
  FOR i=0,n-1 DO BEGIN
    contrib=ld[i].contrib
    k=where(contrib GE max(contrib*trunc_factor))
    IF min(ltemp[k]) LT ltemp_min THEN ltemp_min=min(ltemp[k])
    IF max(ltemp[k]) GT ltemp_max THEN ltemp_max=max(ltemp[k])
  ENDFOR 
ENDELSE   

nt=round((ltemp_max-ltemp_min)/dlogt)+1
ltemp=findgen(nt)*dlogt + ltemp_min

return,ltemp

END
