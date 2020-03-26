EpsilonNeighbourhood:=function(iota)
    local
        M, Bn, dimM, dimBn, image, i, j;

    if not IsHapRegularCWMap(iota) then
        Error("input must be an inclusion of regular CW-complexes");
    fi;

    M:=ShallowCopy(iota!.source!.boundaries);
    Bn:=ShallowCopy(iota!.target!.boundaries);

    dimM:=EvaluateProperty(M,"dimension");
    dimBn:=EvaluateProperty(Bn,"dimension");

    image:=List([1..dimM+1],x->[]);

    for i in [1..dimM+1] do
        for j in [1..Length(M[i])] do
            Add(image[i],iota!.mapping(i-1,j));
        od;
    od;


end;