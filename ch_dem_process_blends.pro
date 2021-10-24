
FUNCTION ch_dem_process_blends, line_data


;+
; NAME:
;     CH_DEM_PROCESS_BLENDS
;
; PURPOSE:
;     Combines any blends flagged in LINE_DATA to create single
;     entries. 
;
; CATEGORY:
;     CHIANTI; differential emission measure.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_PROCESS_BLENDS( Line_Data )
;
; INPUTS:
;     Line_Data: The structure returned by CH_DEM_READ_LINE_IDS.
;
; OUTPUTS:
;     A structure with the same format as LINE_DATA is returned, but
;     with any lines flagged with 'blend_tags' merged into a single
;     feature, with intensities summed. The tags 'index' and 'label'
;     are modified to indicate the composition of the new feature.
;
; EXAMPLE:
;     IDL> line_data=ch_dem_read_line_ids('line_list.txt',int_file='ints.txt')
;     IDL> result=ch_dem_process_blends(line_data)
;
; MODIFICATION HISTORY:
;     Ver.0.1, 31-Jul-2019, Peter Young.
;     Ver.0.2, 9-Aug-2019, Peter Young
;       Added a check to make sure the 'contrib' tag exists. 
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> result = ch_dem_process_blends( line_data )'
  return,-1
ENDIF 

ld=line_data


k=where(trim(ld.blend_tag) NE '',nk)
IF nk NE 0 THEN BEGIN
  tags=ld[k].blend_tag
  tags=tags[uniq(tags,sort(tags))]
  nt=n_elements(tags)
  FOR i=0,nt-1 DO BEGIN
    j=where(ld.blend_tag EQ tags[i],nj)
    IF nj GT 1 THEN BEGIN 
     ;
     ; This is the case where the blend involves lines from the same
     ; ion. Note that I set blend_tag back to an empty string so that
     ; the new blended line is treated as a normal line in the rest of
     ; the code.
     ;
      iref=j[0]
      FOR ii=1,nj-1 DO BEGIN
        ij=j[ii]
        ld[iref].index=ld[iref].index+','+ld[ij].index
        ld[iref].label=ld[iref].label+'+'+ld[ij].label
        ld[iref].int=ld[iref].int+ld[ij].int
        ld[iref].err=sqrt( ld[iref].err^2 + ld[ij].err^2 )
        ld[iref].blend_tag=''
        IF tag_exist(ld,'contrib') THEN ld[iref].contrib=ld[iref].contrib+ld[ij].contrib
        ld[ij].obs_wvl=-1
      ENDFOR 
      ik=where(ld.obs_wvl GT 0)
      ld=ld[ik]
    ENDIF ELSE BEGIN
     ;
     ; Only one line is flagged so I set blend_tag back to an empty
     ; string so that it's processed as a normal line.
     ;
      ld[j].blend_tag=''
    ENDELSE 
  ENDFOR 
ENDIF

return,ld

END
