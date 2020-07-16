################################################################################
############ Input: a CW-subcomplex [X,S] ######################################
################################################################################
########### Output: the corresponding inclusion map f: Y -> X ##################
################################################################################
CWSubcomplexToRegularCWMap:=function(YS)
    local
        map, src, i, j, trg_cell;

    map:={i,j}->YS[2][i+1][j];

    src:=List([1..Length(Filtered(YS[2],y->y<>[]))+1],x->[]);
    src[1]:=List(YS[2][1],x->[1,0]);

    for i in [2..Length(src)-1] do
        for j in [1..Length(YS[2][i])] do
            trg_cell:=YS[1]!.boundaries[i][YS[2][i][j]]*1;
            trg_cell:=trg_cell{[2..trg_cell[1]+1]};
            trg_cell:=List(trg_cell,x->Position(YS[2][i-1],x));
            Add(trg_cell,Length(trg_cell),1);
            Add(src[i],trg_cell*1);
        od;
    od;

    return Objectify(
        HapRegularCWMap,
        rec(
            source:=RegularCWComplex(src),
            target:=ShallowCopy(YS[1]),
            mapping:=map
        )
    );
end;
