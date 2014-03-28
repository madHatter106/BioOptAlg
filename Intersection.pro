; docformat:'rst'

;+
; A procedure to determine the intersection between to sets of unique sorted integers.
;
; :version: 1
;
; :Author:
; 
;   Gert Sclep (European Commission/JRC/IES)
;
; :Copyright:
; 
;   All rights reserved
;   
;-

;+
; :Author:
;   
;   Gert Sclep
; 
; :version: 1
; 
; :History:
;   
;   Created by Gert Sclep: 07/2011
;   
; :Categories:
;   
;   utilities, intersection
;    
; :Description:
; 
;   Obtains the indexes of the intersection between two sorted unique (not containing any doubles) arrays of integers
;   
; :Examples:
; 
;   Example syntax::
;   
;     a = [1,3,5,9,10]
;     b = [0,3,9,12]
;     intersection_exists = intersectionSortedUniqueSets(a,b,index_a,index_b)
;     if (intersection_exists) then begin
;       intersection1 = a[index_a]
;       intersection2 = b[index_b]
;       if (total(intersection1-intersection2) NE 0) then print, 'this can never happen, you are not seeing this'
;       print, 'Intersection:', intersection1
;     endif
;   
; :Pre:
;   
;   `set_a` and `set_b` must be sorted and not contain doubles (i.e. elements must be unique within one set)     
;       
; :Requires:
; 
;   IDL 7.1
;   
; :Returns:
; 
;   0 if both sets have no elements in common
;   
;   1 if both sets have 1 or more elements in common   
;   
; :Params:
; 
;   set_a: in, out, required, type="int/intarr(m)"
;     first array of sorted unique (no doubles) integers
;     
;     if type is scalar, the set is converted to an array of length 1
;   set_b: in, out, required, type="int/intarr(n)"
;     second array of sorted unique (no doubles) integers
;     
;     if type is scalar, the set is converted to an array of length 1   
;   set_a_intersection_indexes: out, type=intarr(o)
;     array of indexes referring to first input array, pointing to the elements of that first array to be found in the intersection
;     
;     if no intersection was found, the output variable will remain undefined
;     
;     the maximum value of the dimension o equals the dimension m of `set_a`
;   set_b_intersection_indexes: out, type=intarr(p)
;     array of indexes referring to second input array, pointing to the elements of that second array to be found in the intersection
;     
;     if no intersection was found, the output variable will remain undefined
;   
;     the maximum value of dimension p equals the dimension n of `set_b` 
;-
FUNCTION intersectionSortedUniqueSets, set_a, set_b, set_a_intersection_indexes, set_b_intersection_indexes
    
    ; If either of the sets is a scalar, make it a vector.
    IF N_ELEMENTS(set_a) EQ 1 && (SIZE(set_a))[0] EQ 0 THEN set_a = [set_a]
    IF N_ELEMENTS(set_b) EQ 1 && (SIZE(set_b))[0] EQ 0 THEN set_b = [set_b]
    
    ; Find the intersection of the ranges.
    mina = set_a[0]
    maxa = set_a[N_ELEMENTS(set_a)-1]
    minb = set_b[0]
    maxb = set_b[N_ELEMENTS(set_b)-1]    
    minab = mina > minb
    maxab = maxa < maxb

    ; If the set ranges don't intersect, leave now.
    IF ((maxa LT minab) AND (minb GT maxab)) OR ((maxb LT minab) AND (mina GT maxab)) THEN BEGIN        
        RETURN, 0
    ENDIF
    
    set_a_restricted_indexes = WHERE(set_a GE minab and set_a LE maxab,nb_set_a_restrict)
    IF nb_set_a_restrict GT 0 THEN BEGIN
      set_a_restricted = set_a[set_a_restricted_indexes]
      set_a_restr_rel = set_a_restricted - minab
    ENDIF ELSE RETURN, 0
    
    set_b_restricted_indexes = WHERE(set_b GE minab and set_b LE maxab,nb_set_b_restrict)
    IF nb_set_b_restrict GT 0 THEN BEGIN
      set_b_restricted = set_b[set_b_restricted_indexes]
      set_b_restr_rel = set_b_restricted - minab
    ENDIF ELSE RETURN, 0
    
    set_a_mapped = LONARR(maxab-minab+1)
    set_a_mapped[set_a_restr_rel] = 1
    
    set_b_mapped = LONARR(maxab-minab+1)
    set_b_mapped[set_b_restr_rel] = 1
    
    sum_set_mapped = set_a_mapped + set_b_mapped
    
    set_a_intersection_indexes = WHERE(sum_set_mapped[set_a_restr_rel] EQ 2,num_intersection_a) + set_a_restricted_indexes[0]
    
    IF num_intersection_a EQ 0 THEN RETURN, 0
    
    set_b_intersection_indexes = WHERE(sum_set_mapped[set_b_restr_rel] EQ 2,num_intersection_b) + set_b_restricted_indexes[0]
    
    RETURN, 1
        
END
