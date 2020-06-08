# inputs a (non-CW) cell complex B, an integer k>=0, and integer lists b, c > 0.
# The function modifies the cell complex B by adding one k-cell e^k whose
# boundary (k-1)-cells are specified by the list b, and whose coboundary
# (k+1)-cells are specified by the list c. Thus the function modifies only the
# lists B[k], B[k+1].
AddCell:=function(B,k,b,c)
end;