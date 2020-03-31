#  Input: an inclusion of regular CW-complexes M -> B^n
# Output: a pair (M*,B^n*) where M* denotes a tubular neighbourhood of M
#         and B^n* denotes a space homotopy equivalent to B^n containing in its
#         interior M*
EpsilonNeighbourhood:=function(iota)
    local
        M, Bn, M_star, Bn_star, dimM,
        dimBn, map, image, i, j, x,
        cocells, k, check;

    if not IsHapRegularCWMap(iota) then
        Error("input must be an inclusion of regular CW-complexes");
    fi;

    M:=ShallowCopy(iota!.source);
    Bn:=ShallowCopy(iota!.target);
    M_star:=List(M!.boundaries,x->[]); # these are
    Bn_star:=Bn!.boundaries*1; # the output spaces

    dimM:=EvaluateProperty(M,"dimension");
    dimBn:=EvaluateProperty(Bn,"dimension");

    map:=List([1..dimM+1],x->[]); # map: M* -> B^n*
    image:=List([1..dimM+1],x->[]); # list of M's cells in Bn's indexing

    for i in [1..dimM+1] do # delete all cells of the original subcomplex
        for j in [1..Length(M!.boundaries[i])] do
            x:=iota!.mapping(i-1,j);
            Add(image[i],x);
            Unbind(Bn_star[i][x]);
        od;
    od;

    for i in [1..dimM+1] do # add as many unique copies of a cell as there
        for j in [1..Length(image[i])] do # were cells in its coboundary
            cocells:=Bn!.coboundaries[i][image[i][j]];
            for k in [2..cocells[1]+1] do
                if i=1 then x:=[1,0]; else x:=[]; fi;
                Add(M_star[i],x);
                Add(Bn_star[i],x);
                Add(map[i],Length(Bn_star[i]));
                
                if IsBound(Bn_star[i+1][cocells[k]]) then
                    check:=Positions(Bn_star[i+1][cocells[k]],image[i][j]);
                    if Length(check)>1 then
                        check:=check[2];
                    else
                        check:=check[1];
                    fi; # reindex the coboundary cells to contain a unique copy
                    Bn_star[i+1][cocells[k]][check]:=Length(Bn_star[i]);
                fi;
            od;
        od;
    od;
    
    return [M_star,Bn_star,map];
end;
