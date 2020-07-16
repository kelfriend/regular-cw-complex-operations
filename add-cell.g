################################################################################
############ Input: a cell complex B, an integer k>=0 and lists of #############
################### positive integers b and c ##################################
################################################################################
########### Output: the cell complex B with an added k-cell whose ##############
################### boundary is specified by b and whose coboundary ############
################### is specified by c ##########################################
################################################################################
AddCell:=function(B,k,b,c)
    local
        i;

    Add(b,Length(b),1);
    Add(B[k+1],b);
    
    for i in [1..Length(c)] do
        Add(B[k+2][c[i]],Length(B[k+1]));
        B[k+2][c[i]][1]:=B[k+2][c[i]][1]+1;
    od;
end);
