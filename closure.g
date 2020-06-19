################################################################################
############ Input: a CW-complex Y and two integers k>=0 & i>1 #################
################################################################################
########### Output: a CW-subcomplex [Y,S] corresponding to the #################
################### topological closure of the ith k-cell of Y #################
################################################################################
ClosureCWCell:=function(Y,k,i)
    local
        clsr, l, m, bnd, n;

    clsr:=List([1..k+1],x->[]);
    Add(clsr[k+1],i);
    
    for l in Reversed([2..k+1]) do
        for m in [1..Length(clsr[l])] do
            bnd:=[];
            for n in [2..Length(Y!.boundaries[l][clsr[l][m]])] do
                Add(bnd,Y!.boundaries[l][clsr[l][m]][n]);
            od;
            clsr[l-1]:=Union(clsr[l-1],bnd);
        od;
    od;

    return [ShallowCopy(Y),clsr];
end;
