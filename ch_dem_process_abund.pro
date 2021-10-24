
FUNCTION ch_dem_process_abund, line_data, abund, ab_elt_fix=ab_elt_fix, $
                               abund_file=abund_file, fixed_abund=fixed_abund, $
                               init=init, out_init=out_init, quiet=quiet, $
                               swtch_ab=swtch_ab

;+
; NAME:
;      CH_DEM_PROCESS_ABUND
;
; CATEGORY:
;      CHIANTI; DEM; abundances.
;
; PURPOSE:
;      Produces a structure that matches the lines in LINE_DATA with
;      abundance information about those lines. The abundance data is
;      used by the CH_DEM suite of routines.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_PROCESS_ABUND( Line_Data )
;
; INPUTS:
;      Line_Data:  A structure in the format returned by
;                  ch_dem_read_line_ids.pro.
;
; OPTIONAL INPUTS:
;      Abund:   A 30-element float array containing the abundances
;               (relative to hydrogen) of the first 30 elements of the
;               periodic table.
;      Ab_Elt_Fix:  Either an integer or a string specifying the
;                   reference element. For example, ab_elt_fix=26 or
;                   'fe' specifies iron. If not specified, then the
;                   routine uses the element with the most number of
;                   lines in LINE_DATA (which is the recommended
;                   option). 
;      Abund_File: An alternative to specifying ABUND is to specify
;                  the name of a CHIANTI format element abundance
;                  file.
;      Init:    When calling ch_dem_process_abund from a DEM
;               routine, then INIT will be the 1D array containing the
;               initial fit parameters of the DEM. For example, if the
;               DEM is a Gaussian (ch_dem_gauss_fit.pro), then init
;               will contain the three initial parameters of the
;               Gaussian.
;      Swtch_Ab:  An array of same size as the number of unique
;                 elements in LINE_DATA. A zero indicates the
;                 element's abundance should be fixed, and a
;                 one indicates the abundance should be a
;                 variable. This is expected to be used only in
;                 special cases. 
;
; KEYWORD PARAMETERS:
;      FIXED_ABUND: If set, then the abundances are not free
;                   parameters. The reference element remains type=1,
;                   but the other parameters become type=3.
;      QUIET:    If set, then information will not be printed to
;                screen. 
;
; OUTPUTS:
;      A structure with N elements, where N is the number of unique
;      elements represented in LINE_DATA. The tags are:
;       .elt_num  The atomic number of the element.
;       .type     Either 0, 1, 2 or 3 (see below).
;       .ind      An IDL list containing the indices of the entries in
;                 LINE_DATA that match the element.
;       .nlines   The number of lines in LINE_DATA that match the
;                 element.
;       .abund    The abundance, relative to hydrogen, of the element.
;       .ratio    The abundance ratio, relative to the reference
;                 element.
;       The "type" tag is used by the DEM software. A type of 1 means
;       the element is the reference element for the DEM. A type of 2
;       means that the element has two or more lines in LINE_DATA and
;       the abundance relative to the reference element will go into
;       the minimization procedure. A type of 0 means there is only
;       one line from the element, and so the relative abundance will
;       not be included in the minimization procedure. A type of 3
;       means the element abundance is fixed relative to the type=1
;       element (/fixed_abund keyword).
;
;       If a problem is found, then a value of -1 is returned.
;
; OPTIONAL OUTPUTS:
;       Out_Init:  If init has been specified, then out_init contains
;                  the new init with the initial abundance parameters
;                  appended to it. For example if Fe is the reference
;                  element, and the Si/Fe ratio is to be varied in the
;                  minimization process then out_init will contain the
;                  initial estimate of the Si/Fe abundance. Note that the
;                  abundances must be specified through abund or
;                  abund_file for this to work. Out_init is used in
;                  the call to the mpfit routines to derived the DEM.
;
; CALLS:
;       READ_ABUND, Z2ELEMENT
;
; EXAMPLES:
;       IDL> line_data=ch_dem_read_line_ids('line_list.txt')
;       IDL> ab=ch_dem_process_abund(line_data)
;       IDL> ab=ch_dem_process_abund(line_data,abund)
;       IDL> ab=ch_dem_process_abund(line_data,ab_elt_fix='fe')
;       IDL> ab=ch_dem_process_abund(line_data,ab_elt_fix=26)
;       IDL> ab=ch_dem_process_abund(line_data,abund_file=!abund_file)
;
;       IDL> init=[1.0,1.0,1.0]
;       IDL> ab=ch_dem_process_abund(line_data,abund_file=!abund_file,
;                       init=init,out_init=out_init) 
;
; MODIFICATION HISTORY:
;       Ver.1, 12-Feb-2019, Peter Young
;       Ver.2, 29-Jul-2019, Peter Young
;          Added /fixed_abund keyword.
;       Ver.3, 16-Dec-2019, Peter Young
;          Added init= and out_init=; updated header; added
;          information text (suppressed with /quiet).
;       Ver.4, 16-Jun-2021, Peter Young
;          Small change to print output.
;       Ver.5, 17-Jun-2021, Peter Young
;          Added SWTCH_AB optional input.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ab=ch_dem_process_abund(line_data [, abund, ab_elt_fix=, abund_file=, /fixed_abund '
  print,'                                      init=, out_init=, /quiet, swtch_ab= ])'
  return,-1
ENDIF 
  


elt_num=line_data.element_num
elt_num_unq=elt_num[uniq(elt_num,sort(elt_num))]
n_elt=n_elements(elt_num_unq)

;
; Do some checks on swtch_ab.
;
ns=n_elements(swtch_ab)
IF ns NE 0 THEN BEGIN
   IF ns NE n_elt THEN BEGIN
      print,'% CH_DEM_PROCESS_ABUND: the input SWTCH_AB should have the same number of elements as there are '
      print,'                        unique elements in the line list. Please check your inputs. Returning...'
      print,'   No. of elements (swtch_ab):  '+trim(ns)
      print,'   No. of elements (line list): '+trim(n_elt)
      return,-1
   ENDIF
  ;
   k=where(swtch_ab EQ 0,nk)
   IF nk EQ 0 THEN BEGIN
      print,"% CH_DEM_PROCESS_ABUND: the input SWTCH_AB must contain at least one zero value (which specifies"
      print,"                        that an element's abundance is fixed). Returning..."
      return,-1
   ENDIF 
ENDIF 

;
; Create the output structure and populate some of the tags.
;
str={elt_num: 0,$
     type: 0, $
     init_ab_ind: -1, $
     abund: -1., $
     error: -1., $
     ratio: -1., $
     ratio_error: -1.}
abstr=replicate(str,n_elt)
abstr.elt_num=elt_num_unq
nlines=intarr(n_elt)
FOR i=0,n_elt-1 DO BEGIN
  k=where(abstr[i].elt_num EQ line_data.element_num,nk)
  nlines[i]=nk
ENDFOR 



;
; Use AB_ELT_FIX to specify the reference element (type=1), or
; automatically choose it.
;
IF n_elements(ab_elt_fix) NE 0 THEN BEGIN
   z2element,indgen(30)+1,elt,/symbol
   IF datatype(ab_elt_fix) EQ 'STR' THEN BEGIN
      chck=where(strlowcase(ab_elt_fix) EQ strlowcase(elt),nchck)
      elt_num=chck[0]+1
      k=where(elt_num EQ abstr.elt_num,nk)
      IF nk EQ 0 THEN BEGIN
         print,'% CH_DEM_PROCESS_ABUND: the element name contained in AB_ELT_FIX does not match any of the lines in LINE_DATA. Returning...'
         return,-1
      ENDIF
   ENDIF ELSE BEGIN
      k=where(abstr.elt_num EQ ab_elt_fix,nk)
      IF nk EQ 0 THEN BEGIN
         print,'% CH_DEM_PROCESS_ABUND: the element number contained in AB_ELT_FIX does not match any of the lines in LINE_DATA. Returning...'
         return,-1
      ENDIF
   ENDELSE
  ;
  ; Need to check there isn't a conflict with SWTCH_AB.
  ;
   IF ns GT 0 THEN BEGIN
      IF swtch_ab[k[0]] EQ 1 THEN BEGIN
         print,'% CH_DEM_PROCESS_ABUND: AB_ELT_FIX specifies an element to be fixed, but SWTCH_AB specifies the'
         print,'                        same element to be variable. I suggest not specifying AB_ELT_FIX. '
         print,'                        Returning...'
         return,-1
      ENDIF 
   ENDIF 
   abstr[k[0]].type=1
ENDIF ELSE BEGIN
  ;
  ; If swtch_ab has been specified, then I need to choose the
  ; reference element from those with swtch_ab=0
  ;
   IF ns EQ 0 THEN BEGIN 
      getmax=max(nlines,imax)
      abstr[imax].type=1
   ENDIF ELSE BEGIN
      k=where(swtch_ab EQ 0)
      nlines2=nlines[k]
      getmax=max(nlines2,imax)
      abstr[k[imax]].type=1
   ENDELSE 
ENDELSE

;
; Identify any type=2 elements. The abundances for these are free to
; vary in the minimization procedure.
;
k=where(nlines GT 1 AND abstr.type NE 1,nk)
IF nk NE 0 THEN abstr[k].type=2


;
; The abundances can be specified directly (through ABUND) or through
; an abundance filename (ABUND_FILE). If neither, then ABUND remains
; undefined. 
;
IF n_elements(abund) EQ 0 AND n_elements(abund_file) NE 0 THEN BEGIN
  chck=file_info(abund_file)
  IF chck.exists THEN read_abund,abund_file,abund,ref
ENDIF
;
IF n_elements(abund) NE 0 THEN BEGIN
  abstr.abund=abund[abstr.elt_num-1]
  k=where(abstr.type EQ 1)
  abstr.ratio=abstr.abund/abstr[k].abund
ENDIF 


;
; Implementation of swtch_ab
;
IF ns GT 0 THEN BEGIN
  ;
  ; This switches any type 0 or type 2 elements to type 3 if
  ; swtch_ab=0. That is, instead of being treated as free parameters,
  ; they will be fixed.
  ;
   k=where(abstr.type EQ 0 OR abstr.type EQ 2,nk)
   IF nk NE 0 THEN BEGIN
      j=where(swtch_ab[k] EQ 0,nj)
      IF nj NE 0 THEN abstr[k[j]].type=3
   ENDIF 
ENDIF 

;
; This keyword makes the abundances fixed in the minimization
; process. There is still a reference element (type=1), but all other
; elements are considered type=3. See ch_dem_gauss_fit_fn for how type
; is interpreted by the fit routine.
;
; Note that this over-rides SWTCH_AB!
;
IF keyword_set(fixed_abund) THEN BEGIN
  k=where(abstr.type NE 1,nk)
  IF nk NE 0 THEN abstr[k].type=3
ENDIF 

;
; line_data.ab_ind is the abstr index of the element in line_data.
;
FOR i=0,n_elt-1 DO BEGIN
  k=where(line_data.element_num EQ elt_num_unq[i])
  line_data[k].ab_ind=i
ENDFOR 



;
; If INIT has been specified and the abundance array is defined, then
; modify INIT to include the abundances that will be varied.
;
IF n_elements(abund) NE 0 AND n_elements(init) NE 0 THEN BEGIN
  ninit=n_elements(init)
 ;
 ; Only "type 2" elements are added to the initial parameter array. If
 ; none, then just return the original array.
 ;
  k=where(abstr.type EQ 2,nk)
  IF nk NE 0 THEN BEGIN 
    FOR i=0,nk-1 DO BEGIN
      j=k[i]
      out_init=[init,abstr[j].ratio]
      ni=n_elements(out_init)
      abstr[j].init_ab_ind=ni-1
    ENDFOR
  ENDIF ELSE BEGIN
    out_init=init
  ENDELSE 
ENDIF

IF NOT keyword_set(quiet) THEN BEGIN
   n_ab=n_elements(abstr)
   print,'% CH_DEM_PROCESS_ABUND: there are '+trim(n_ab)+' different elements. Their types are:'
   FOR i=0,n_ab-1 DO BEGIN
      type=abstr[i].type
      CASE type OF
         0: type_text='derived algebraically (single line)'
         1: type_text='reference'
         2: type_text='variable (multiple lines)'
         3: type_text='fixed (/fixed_abund keyword)'
      ENDCASE
      type_text=strpad(type_text,40,fill=' ',/after)
      z2element,abstr[i].elt_num,name,/symbol
      name=strpad(name,2,fill=' ',/after)
      print,format='(5x,i3,".",i5,2x,a2,"  type=",i1,"  -- ",a40)',i+1,abstr[i].elt_num,name,abstr[i].type,type_text
   ENDFOR
   print,''
ENDIF 

return,abstr

END
