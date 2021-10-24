
FUNCTION ch_wgfa_select, ion_name, wmin=wmin, wmax=wmax, all=all

;+
; NAME:
;     CH_WGFA_SELECT
;
; PURPOSE:
;     Allows a user to select lines from the list within the .wgfa
;     file, and returns the wgfa structure containing only these
;     lines. 
;
; CATEGORY:
;     CHIANTI; utilities.
;
; CALLING SEQUENCE:
;     Result = CH_WGFA_SELECT( Ion_Name )
;
; INPUTS:
;     Ion_Name: The name of an ion in CHIANTI format (e.g., 'o_6')
;
; OPTIONAL INPUTS:
;     Wmin:  Show only wavelengths > WMIN (angstroms).
;     Wmax:  Show only wavelengths < WMAX (angstroms).
;	
; KEYWORD PARAMETERS:
;     ALL:   By default only lines with observed wavelengths are shown
;            in the widget. This keyword allows lines with theoretical
;            wavelengths to be shown as well.
;
; OUTPUTS:
;     Returns a structure in the same format as that returned by
;     READ_WGFA_STR, but containing only those lines selected by the
;     user.
;
;     If no lines are selected, then -1 is returned.
;
; CALLS:
;     CONVERTNAME, ZION2FILENAME, READ_ELVLC, READ_WGFA_STR,
;     CH_XMENU_SEL. 
;
; EXAMPLE:
;     IDL> wgfa=ch_wgfa_select('o_6')
;     IDL> wgfa=ch_wgfa_select('o_6',wmin=1030,wmax=1040)
;     IDL> wgfa=ch_wgfa_select('o_6',wmin=1030,wmax=1040,/all)
;
; MODIFICATION HISTORY:
;     Ver.1, 9-Apr-2019, Peter Young
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> wgfa=ch_wgfa_select( ion_name [, wmin=, wmax=, /all] )'
  return,-1
ENDIF 

convertname,ion_name,iz,ion
zion2filename,iz,ion,fname

read_wgfa_str,fname+'.wgfa',wgfa
read_elvlc,fname+'.elvlc',elvlc=elvlc

IF n_elements(wmin) EQ 0 THEN wmin=min(abs(wgfa.wvl))
IF n_elements(wmax) EQ 0 THEN wmax=max(abs(wgfa.wvl))
k=where(abs(wgfa.wvl) GE wmin AND abs(wgfa.wvl) LE wmax,nk)
IF nk EQ 0 THEN BEGIN
  print,'% CH_WGFA_SELECT: there are no lines between WMIN and WMAX. Returning...'
  return,-1
ENDIF
wgfa=wgfa[k]
i=sort(abs(wgfa.wvl))
wgfa=wgfa[i]

IF NOT keyword_set(all) THEN BEGIN
  k=where(wgfa.wvl GT 0.,nk)
  IF nk EQ 0 THEN BEGIN
    print,'% CH_WGFA_SELECT: there are no lines between WMIN and WMAX. Returning...'
    return,-1
  ENDIF
  wgfa=wgfa[k]
ENDIF 

;
; Create the text list that the user will see in the widget.
;
options=strarr(nk)
FOR i=0,nk-1 DO BEGIN
  l1=wgfa[i].lvl1
  l2=wgfa[i].lvl2
  desc1=elvlc.data[l1-1].full_level
  desc2=elvlc.data[l2-1].full_level
 ;
  len=strlen(trim(string(ceil(max(abs(wgfa.wvl))))))
  wformat='f'+trim(len+4)+'.3'
 ;
  len=strlen(trim(max(wgfa.lvl2)))
  lformat='i'+trim(len+1)
 ;
  IF wgfa[i].wvl LT 0 THEN text_add=' *' ELSE text_add='  '
  options[i]=string(format='('+wformat+',1x,a1,a2,1x,'+lformat+',1x,a1,1x,'+lformat+')', $
                    abs(wgfa[i].wvl),string(197b),text_add,l1,'-',l2)
  options[i]=options[i]+'  '+desc1+' - '+desc2
ENDFOR 

index = ch_xmenu_sel(options, tit=' Select lines ', $
                     text=['If more than one line is selected, ',$
                           ' the emissivities of the lines will be summed.', $
                           'A  * indicates that the wavelength is theoretical' ])

IF index[0] EQ -1 THEN return,index[0] ELSE return,wgfa[index]

return,wgfa

END
