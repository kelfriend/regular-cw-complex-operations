################################################################################
############ Input: a strictly cellular map f: Y -> X of regular ###############
################### CW-complexes ###############################################
################################################################################
########### Output: the CW-subcomplex [X,S] corresponding to the ###############
################### image of f #################################################
################################################################################
RegularCWMapToCWSubcomplex:=function(f)
    local
        dim, src_bnd, S, i, j;

    dim:=EvaluateProperty(f!.source,"dimension")*1;
    src_bnd:=f!.source!.boundaries*1;
    S:=List([0..dim],x->[]);

    for i in [1..dim+1] do
        for j in [1..Length(src_bnd[i])] do
            Add(S[i],f!.mapping(i-1,j));
        od;
    od;
    
    return [ShallowCopy(f!.target),S];
end;
