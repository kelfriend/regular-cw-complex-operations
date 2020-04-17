#  Input: an inclusion of regular CW-complexes M -> B^n
# Output: a pair (M*,B^n*) where M* denotes a tubular neighbourhood of M
#         and B^n* denotes a space homotopy equivalent to B^n containing in its
#         interior M*
SubspaceComplement:=function(iota)
    local
        M, Bn, M_star, Bn_star, dimM,
        dimBn, map, image, i, j, x, replicate,
        cocells, k, check, cell, gap, original,
        bndbnd, l, int, m, handicap, Bn_final,
        count, bool;

    if not IsHapRegularCWMap(iota) then
        Error("input must be an inclusion of regular CW-complexes");
    fi;

    M:=ShallowCopy(iota!.source);
    Bn:=ShallowCopy(iota!.target);
    M_star:=List(M!.boundaries*1,x->[]); # these are the output spaces
    Bn_star:=List(Bn!.boundaries*1,x->Concatenation(x,[[0]]));
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
    replicate:=List([1..dimM+1],x->List([1..Length(Bn!.boundaries[x])],y->[]));

    for i in [1..dimM+1] do # add as many unique copies of a cell as there
        for j in [1..Length(image[i])] do # were cells (not in the subcomplex)
            cocells:=Bn!.coboundaries[i][image[i][j]]*1; # in its coboundary
            cocells:=cocells{[2..Length(cocells)]};
            cocells:=Filtered(cocells,x->IsBound(Bn_star[i+1][x]));
            for k in [1..Length(cocells)] do
                if i=1 then
                    x:=[1,0];
                else
                    x:=[];
                fi;
                Add(M_star[i],x);
                Add(Bn_star[i],x);
                Add(map[i],Length(Bn_star[i]));
                Add(replicate[i][image[i][j]],Length(Bn_star[i]));
                if IsBound(Bn_star[i+1][cocells[k]]) and
                   Bn_star[i+1][cocells[k]]<>[0] then
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

    for i in [3..dimBn+1] do # detect cells with empty boundary and
        for j in [1..Length(Bn_star[i])] do # attempt to patch them
            if IsBound(Bn_star[i][j]) then
                cell:=Bn_star[i][j];
                if cell<>[0] then
                    bool:=true;
                fi;
            else
                bool:=false;
            fi;
            gap:=false;
            bndbnd:=[];
            if bool then
                for k in cell{[2..Length(cell)]} do
                    if Bn_star[i-1][k]=[] then
                        gap:=true;
                    fi;
                    for l in Bn_star[i-1][k]{[2..Length(Bn_star[i-1][k])]} do
                        Add(bndbnd,l);
                    od;
                od;
                bndbnd:=Filtered(bndbnd,x->Length(Positions(bndbnd,x))=1);
                if gap then
                    for k in Filtered(cell{[2..Length(cell)]},x->Bn_star[i-1][x]=[])
                        do # find out where the empty (i-1)-cell was rep'd from
                        original:=Position(List(replicate[i-1]*1,x->k in x),true);
                        original:=Bn!.boundaries[i-1][original];
                        original:=original{[2..Length(original)]};
                        for l in original do
                            int:=Intersection(replicate[i-2][l],bndbnd);
                            if int<>[] then
                                for m in [1..Length(int)] do
                                    Add(Bn_star[i-1][k],int[m]);
                                od;
                            fi;
                        od;
                        Add(Bn_star[i-1][k],Length(Bn_star[i-1][k]),1);
                    od;
                fi;
            fi;
        od;
    od;
    for i in [3..dimBn+1] do # detect gaps in the boundary of the boundary in
        for j in [1..Length(Bn_star[i])] do # the chain complex and patch them
            if IsBound(Bn_star[i][j]) then
                cell:=Bn_star[i][j];
                if cell<>[0] then
                    bool:=true;
                fi;
            else
                bool:=false;
            fi;
            bndbnd:=[];
            if bool then
                for k in cell{[2..Length(cell)]} do
                    for l in Bn_star[i-1][k]{[2..Length(Bn_star[i-1][k])]} do
                        Add(bndbnd,l);
                    od;
                od;
                bndbnd:=Filtered(bndbnd,x->Length(Positions(bndbnd,x))=1);
                if bndbnd<>[] then
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
            fi;
        od;
    od;

    handicap:=List(Bn_star*1,x->[]); # reindex everything so that Bn_final
    Bn_final:=handicap*1; # yields a regular CW-complex from Bn_star
    for i in [1..Length(Bn_star)] do
        for j in [1..Length(Bn_star[i])] do
            if IsBound(Bn_star[i][j]) then
                if Bn_star[i][j]=[0] then
                    Unbind(Bn_star[i][j]);
                fi;
            fi;
        od;
    od;

    for i in [1..Length(Bn_star)] do
        for j in [1..Length(Bn_star[i])] do
            if IsBound(Bn_star[i][j]) then
                if Bn_star[i][j]<>[0] then
                    if i=1 then
                        Add(Bn_final[1],[1,0]);
                    else
                        cell:=Bn_star[i][j]*1;
                        if cell<>[0] then
                            for k in [2..Length(cell)] do
                                count:=0;
                                for l in [1..cell[k]] do
                                    if not IsBound(Bn_star[i-1][l]) then
                                        count:=count+1;
                                    fi;
                                od;
                                cell[k]:=cell[k]-count;
                            od;
                            Add(Bn_final[i],cell);
                        fi;
                    fi;
                fi;
            fi;
        od;
    od;

    return Bn_star;#Bn_final;#[M_star,Bn_star,map];
        #Objectify(
        #    HapRegularCWMap,
        #    rec(
        #        source:=M_star,
        #        target:=Bn_final,
        #        mapping:=map
        #    )
        #);
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
m2:=[[[1,0],[1,0]],[[2,1,2]],[]];
b22:=[
    List([1..6],x->[1,0]),
    [ [2,1,2], [2,1,3], [2,1,4], [2,2,5], [2,1,6], [2,3,4], [2,4,5], [2,5,6], [2,3,6] ],
    [ [3,2,5,9], [3,2,3,6], [4,1,4,5,8], [4,1,4,3,7] ],
    []
];
mp2:=function(i,j)
if i=1 then return 1; fi;
if i=0 then return j; fi;
end;
omicron:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(m2),
        target:=RegularCWComplex(b22),
        mapping:=mp2
    )
);
m3:=[[[1,0],[1,0]],[[2,1,2]],[]];
b23:=[
    List([1..7],x->[1,0]),
    [[2,1,2],[2,1,2],[2,1,5],[2,1,6],[2,2,3],[2,2,6],
    [2,3,4],[2,3,4],[2,3,6],[2,4,5],[2,4,7],[2,5,6],
    [2,5,7],[2,6,7]],
    [[2,1,2],[3,2,4,6],[3,5,6,9],[2,7,8],[4,8,9,11,14],
    [3,10,11,13],[3,12,13,14],[3,3,4,12]],
    []
];
mp3:=function(i,j)
if i=0 then if j=1 then return 6; else return 7; fi; else return 14; fi;
end;
phi:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(m3),
        target:=RegularCWComplex(b23),
        mapping:=mp3
    )
);
m4:=[[[1,0],[1,0]],[[2,1,2]],[]];
b3:=[
    List([1..8],x->[1,0]),
    [[2,1,2],[2,2,3],[2,3,4],[2,4,5],[2,5,6],[2,1,6],[2,1,7],[2,7,8],[2,4,8],
    [2,2,7],[2,6,7],[2,3,8],[2,5,8]],
    [[3,1,7,10],[3,6,7,11],[4,11,5,8,13],[4,2,10,12,8],[3,3,12,9],[3,9,13,4],
    [6,1,2,3,4,5,6],[6,1,2,3,4,5,6]],
    [[7,1,2,3,4,5,6,7],[7,1,2,3,4,5,6,8]],
    []];
mp4:=function(i,j)
    if i=0 then
        if j=1 then
            return 7;
        else
            return 8;
        fi;
    else
        return 8;
    fi;
end;
psi:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(m4),
        target:=RegularCWComplex(b3),
        mapping:=mp4
    )
);
