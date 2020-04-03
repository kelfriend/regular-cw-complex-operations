#  Input: an inclusion of regular CW-complexes M -> B^n
# Output: a pair (M*,B^n*) where M* denotes a tubular neighbourhood of M
#         and B^n* denotes a space homotopy equivalent to B^n containing in its
#         interior M*
EpsilonNeighbourhood:=function(iota)
    local
        M, Bn, M_star, Bn_star, dimM,
        dimBn, map, image, i, j, x, replicate,
        cocells, k, check, cell, gap1,
        gap2, bndbnd, l;

    if not IsHapRegularCWMap(iota) then
        Error("input must be an inclusion of regular CW-complexes");
    fi;

    M:=ShallowCopy(iota!.source);
    Bn:=ShallowCopy(iota!.target);
    M_star:=List(M!.boundaries,x->[]); # these are
    Bn_star:=Bn!.boundaries*1; # the output spaces

    dimM:=EvaluateProperty(M,"dimension");
    dimBn:=EvaluateProperty(Bn,"dimension");

    map:=List([1..dimM+2],x->[]); # map: M* -> B^n*
    image:=List([1..dimM+1],x->[]); # list of M's cells in Bn's indexing

    for i in [1..dimM+1] do # delete all cells of the original subcomplex
        for j in [1..Length(M!.boundaries[i])] do
            x:=iota!.mapping(i-1,j);
            Add(image[i],x);
            Unbind(Bn_star[i][x]);
        od;
    od;
    replicate:=List([1..dimM+1],x->List([1..Length(Bn_star[x])],y->[]));

    for i in [1..dimM+1] do # add as many unique copies of a cell as there
        for j in [1..Length(image[i])] do # were cells in its coboundary
            cocells:=Bn!.coboundaries[i][image[i][j]];
            for k in [2..cocells[1]+1] do
                if i=1 then x:=[1,0]; else x:=[]; fi;
                Add(M_star[i],x);
                Add(Bn_star[i],x);
                Add(map[i],Length(Bn_star[i]));
                Add(replicate[i][image[i][j]],Length(Bn_star));
                
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

    for i in [3..dimBn+1] do
        for j in [1..Length(Bn_star[i])] do
            cell:=Bn_star[i][j];
            gap1:=false;
            gap2:=false;
            bndbnd:=[];
            for k in cell{[2..Length(cell)]} do
                if Bn_star[i-1][k]=[] then
                    gap1:=true;
                fi;
                for l in Bn_star[i-1][k]{[2..Length(Bn_star[i-1][k])]} do
                    Add(bndbnd,l);
                od;
            od;
            bndbnd:=Filtered(bndbnd,x->Length(Positions(bndbnd,x))=1);
            if bndbnd<>[] and not gap1 then
                gap2:=true;
            fi;
            if gap1 then
                for k in Filtered(cell{[2..Length(cell)]},x->Bn_star[i-1][x]=[])
                    do # find out where the empty (n-1)-cell was rep'd from
                    original:=Position(List(replicate[i-1],x->k in x),true);
                    # now find the duplicates of that original cells boundary
                od;
            elif gap2 then
                bndbnd:=Concatenation([Length(bndbnd)],SortedList(bndbnd));
                
                Add(Bn_star[i-1],bndbnd);
                Add(M_star[i-1],"*");
                Add(map[i-1],Length(Bn_star[i-1]));

                Add(cell,Length(Bn_star[i-1]));
                Bn_star[i][j]:=Concatenation(
                    [Length(cell)-1],
                    cell{[2..Length(cell)]}
                );
            fi;
        od;
    od;    
    
    return [M_star,Bn_star,map];
end;
m:=[ [ [1,0] ], [ ] ];
b2:=[
    [ [1,0], [1,0], [1,0], [1,0], [1,0] ],
    [ [2,1,2], [2,1,3], [2,1,4], [2,1,5], [2,2,3], [2,3,4], [2,4,5], [2,2,5] ],
    [ [3,1,4,8], [3,1,2,5], [3,2,3,6], [3,3,4,7] ],
    [ ]
];
mp:={i,j}->1;
iota:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(m),
        target:=RegularCWComplex(b2),
        mapping:=mp
    )
);
