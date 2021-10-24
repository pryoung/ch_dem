
FUNCTION ch_dem_mcmc_emis, line_data, ltemp, ldens=ldens, lpress=lpress

;+
; NAME:
;	ROUTINE_NAME
;
; PURPOSE:
;     Computes the CHIANTI 'emis' array needed by mcmc_dem routine.
;
; CATEGORY:
;	Put a category (or categories) here.  For example:
;	Widgets.
;
; CALLING SEQUENCE:
;	Write the calling sequence here. Include only positional parameters
;	(i.e., NO KEYWORDS). For procedures, use the form:
;
;	ROUTINE_NAME, Parameter1, Parameter2, Foobar
;
;	Note that the routine name is ALL CAPS and arguments have Initial
;	Caps.  For functions, use the form:
; 
;	Result = FUNCTION_NAME(Parameter1, Parameter2, Foobar)
;
; INPUTS:
;	Parm1:	Describe the positional input parameters here. Note again
;		that positional parameters are shown with Initial Caps.
;
; OPTIONAL INPUTS:
;	Parm2:	Describe optional inputs here. If you don't have any, just
;		delete this section.
;	
; KEYWORD PARAMETERS:
;	KEY1:	Document keyword parameters like this. Note that the keyword
;		is shown in ALL CAPS!
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;	Describe optional outputs here.  If the routine doesn't have any, 
;	just delete this section.
;
; COMMON BLOCKS:
;	BLOCK1:	Describe any common blocks here. If there are no COMMON
;		blocks, just delete this entry.
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; CALLS:
;     CH_DEM_PROCESS_INDEX_STRING, CH_LOOKUP_EMISS
;
; EXAMPLE:
;	Please provide a simple example here. An example from the
;	DIALOG_PICKFILE documentation is shown below. Please try to
;	include examples that do not rely on variables or data files
;	that are not defined in the example code. Your example should
;	execute properly if typed in at the IDL command line with no
;	other preparation. 
;
; MODIFICATION HISTORY:
; 	Written by:	Your name here, Date.
;	July, 1994	Any additional mods get described here.  Remember to
;			change the stuff above if you add a new keyword or
;			something!
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> emis = ch_dem_mcmc_emis( line_data, ltemp [, ldens=, lpress= ] )'
  return,-1
ENDIF 

IF n_elements(line_data) GT 1 THEN BEGIN
  print,'% CH_DEM_MCMC_EMIS: The input LINE_DATA must be a single-element structure. Returning...'
  return,-1
ENDIF

indstr=ch_dem_process_index_string(line_data.index)
n=n_elements(indstr)

nt=n_elements(ltemp)
emis=dblarr(nt)

IF n_elements(lpress) NE 0 THEN pressure=10.^lpress

em=ch_lookup_emiss(line_data.ion,ltemp=ltemp,ldens=ldens,pressure=pressure,/quiet)

FOR i=0,n-1 DO BEGIN
  k=where(em.level1 EQ indstr[i].lower AND em.level2 EQ indstr[i].upper,nk)
  IF nk NE 0 THEN BEGIN
    emis=emis+em[k[0]].em*1d23
  ENDIF 
ENDFOR

junk=temporary(em)

return,emis

END
