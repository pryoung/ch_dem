

FUNCTION ch_dem_compute_ints, data


;
; Unpack the data structure
;
phi=data.phi
line_data=data.line_data
ab_ref=data.ab_ref
abstr=data.abstr
ltemp=data.ltemp
p=data.p

dlogt=ltemp[1]-ltemp[0]
t=10.^ltemp


nl=n_elements(line_data)
int=dblarr(nl)

;
; Now compute the intensities. Note how abfac is computed differently
; depending on the abundance type (ab_type).
; 
;; FOR i=0,nl-1 DO BEGIN
;;   CASE ab_type[i] OF
;;     1: abfac=1.0 
;;     2: abfac=p[init_ab_ind[i]]
;;     3: abfac=abund[i]/ab_ref
;;   ENDCASE 
;;   int[i]=total(phi * abfac * ab_ref * contrib_data[*,i] * dlogt * t* alog(10.) )
;; ENDFOR 

ab_type=abstr[line_data.ab_ind].type
init_ab_ind=abstr[line_data.ab_ind].init_ab_ind
abund=abstr[line_data.ab_ind].abund
FOR i=0,nl-1 DO BEGIN
  CASE ab_type[i] OF
    1: abfac=1.0 
    2: abfac=p[init_ab_ind[i]]
    3: abfac=abund[i]/ab_ref
  ENDCASE
  contrib=line_data[i].contrib
  int[i]=total(phi * abfac * ab_ref * contrib * dlogt * t* alog(10.) )
ENDFOR 


;
; New
;
;; FOR i=0,nl-1 DO BEGIN
;;   nlines=line_data[i].nlines
;;   FOR j=0,nlines-1 DO BEGIN
;;     elt_num=line_data[i].element_num[j]
;;     k=where(abstr.elt_num EQ elt_num)
;;     ab_type=abstr[k].type
;;     abund=abstr[k].abund
;;     init_ab_ind=abstr[k].init_ab_ind
;;     CASE ab_type OF
;;       1: abfac=1.0 
;;       2: abfac=p[init_ab_ind]
;;       3: abfac=abund/ab_ref
;;     ENDCASE
;;     contrib=line_data[i].contrib[j]
;;     int[i]=int[i] + total(phi * abfac * ab_ref * contrib * dlogt * t* alog(10.) )
;;   ENDFOR 
;; ENDFOR 

;; print,format='(10f8.2)',p

;
; Delete data since it's not needed any more
;
data=0.

return,int

END
