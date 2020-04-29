SubspaceComplement:=function(inc)
    local
        src, trg, dim_s, dim_t, SRC, TRG,
        map, i, j, rep, cell, cobnd, k, x,
        cocell, pos, face, bndbnd, l, base,
        int, len, edge, bnd, original, pairs,
        match, TRG_final, map_final, count;

    src:=ShallowCopy(inc!.source);
    trg:=ShallowCopy(inc!.target);

    dim_s:=EvaluateProperty(src,"dimension");
    dim_t:=EvaluateProperty(trg,"dimension");

    SRC:=List(src!.boundaries*1,x->[]);
    TRG:=List(trg!.boundaries*1,x->Concatenation(x,[[0]]));

    map:=List([1..dim_s+1],x->[]);

    # step 1: delete all cells of src from trg
    for i in [0..dim_s] do
        for j in [1..Length(src!.boundaries[i+1])] do
            Unbind(TRG[i+1][inc!.mapping(i,j)]);
        od;
    od;

    # step 2: replicate each deleted src cell as many times as there were
    # cells of trg \ src in its coboundary
    rep:=List([1..dim_s+1],x->List([1..Length(trg!.boundaries[x])],y->[]));
    for i in [0..dim_s] do
        for j in [1..Length(src!.boundaries[i+1])] do
            cell:=inc!.mapping(i,j);

            cobnd:=trg!.coboundaries[i+1][cell];
            cobnd:=cobnd{[2..cobnd[1]+1]};
            cobnd:=Filtered(cobnd,x->IsBound(TRG[i+2][x]));

            for k in [1..Length(cobnd)] do
                if i=0 then
                    x:=[1,0];
                else
                    x:=[];
                fi;
                Add(SRC[i+1],x);
                Add(TRG[i+1],x);

                Add(map[i+1],Length(TRG[i+1]));
                Add(rep[i+1][inc!.mapping(i,j)],Length(TRG[i+1]));

                # now to reindex the boundaries of the (n+1)-cells in the
                # coboundary so that they each contain a unique copy of the
                # deleted n-cell
                if IsBound(TRG[i+2][cobnd[k]]) then
                    cocell:=TRG[i+2][cobnd[k]];
                    pos:=Positions(cocell,cell)[Length(Positions(cocell,cell))];
                    cocell[pos]:=Length(TRG[i+1]);
                fi;
            od;
        od;
    od;

    # step 3: patch all `gaps of type 1', i.e., determine the boundary
    # of those cells which were added as replicated cells
    for i in [3..dim_t+1] do
        for j in [1..Length(TRG[i])] do
            if IsBound(TRG[i][j]) then
                if TRG[i][j]<>[] then
                    face:=TRG[i][j]*1;
                    face:=face{[2..face[1]+1]};

                    bndbnd:=List(face,x->TRG[i-1][x]);
                    for k in [1..Length(bndbnd)] do
                        if bndbnd[k]<>[] then
                            bndbnd[k]:=bndbnd[k]{[2..bndbnd[k][1]+1]};
                        fi;
                    od;
                    bndbnd:=Concatenation(bndbnd);
                    for k in face do
                        if TRG[i-1][k]=[] then
                            for l in [1..Length(rep[i-1])] do
                                if k in rep[i-1][l] then
                                    base:=trg!.boundaries[i-1][l]*1;
                                    base:=base{[2..Length(base)]};
                                    base:=List(base,x->rep[i-2][x]);
                                    base:=Concatenation(base);

                                    int:=Intersection(bndbnd,base);
                                    if i=3 then
                                        len:=Length(int);
                                        if len=1 then
                                            Add(TRG[1],[1,0]);
                                            Add(int,Length(TRG[1]));
                                            Add(int,2,1);
                                            TRG[2][k]:=int*1;
                                        elif len=0 then
                                            Add(TRG[1],[1,0]);
                                            Add(TRG[1],[1,0]);
                                            TRG[2][k]:=[
                                                2,
                                                Length(TRG[1])-1,
                                                Length(TRG[1])
                                            ];
                                        else
                                            Add(int,2,1);
                                            TRG[2][k]:=int*1;
                                        fi;
                                    else
                                        Add(int,Length(int),1);
                                        TRG[i-1][k]:=int*1;
                                    fi;
                                fi;
                            od; 
                        fi;
                    od;
                fi;
            fi;
        od;
    od;

    # step 4: patch all `gaps of type 2', i.e., add cells where necessary in
    # order to ensure that the boundary of the boundary is 0 in the chain
    # complex
    for i in [3..dim_t+1] do
        for j in [1..Length(TRG[i])*1] do
            if IsBound(TRG[i][j]) then
                if TRG[i][j]<>[] then
                    bndbnd:=TRG[i][j]*1;
                    bndbnd:=bndbnd{[2..bndbnd[1]+1]};
                    bndbnd:=List(bndbnd*1,x->TRG[i-1][x]{[2..TRG[i-1][x][1]+1]});
                    bnd:=bndbnd*1;
                    bndbnd:=Concatenation(bndbnd);
                    bndbnd:=Filtered(bndbnd,x->Length(Positions(bndbnd,x))<2);
                    if bndbnd<>[] then
                        if i=3 and Length(bndbnd)>2 then
                            
                            #original:=[];
                            #for k in [2..TRG[3][j][1]+1] do
                            #    for l in [1..Length(rep[2])] do
                            #        if TRG[3][j][k] in rep[2][l] then
                            #            Add(original,l);
                            #        fi;
                            #    od;
                            #od;
                            #pairs:=[];
                            #for k in original do
                            #    for l in original do
                            #        if k<>l then
                            #            int:=Intersection(
                            #                trg!.boundaries[2][k]{[2,3]},
                            #                trg!.boundaries[2][l]{[2,3]}
                            #            );
                            #            if int<>[] then
                            #                Add(
                            #                    pairs,
                            #                    [rep[2][k],rep[2][l]]
                            #                );
                            #            fi;
                            #        fi;
                            #    od;
                            #od;
                            #pairs:=Set(List(pairs,Set));
                            #for k in [1..Length(pairs)] do
                            #    match:=List(
                            #        pairs[k]*1,
                            #        x->List(x,y->TRG[2][y]{[2,3]})
                            #    );
                            #    match:=List(
                            #        match,
                            #        x->List(
                            #            x,
                            #            y->Filtered(
                            #                y,
                            #                z->z in bndbnd
                            #            )
                            #        )
                            #    );
                            #od;
                        else
                            Add(bndbnd,Length(bndbnd),1);
                            Add(TRG[i-1],bndbnd);

                            Add(TRG[i][j],Length(TRG[i-1]));
                            TRG[i][j][1]:=TRG[i][j][1]+1;
                        fi;
                    fi;
                fi;
            fi;
        od;
    od;
if false then
    # step 5: reindex TRG & map so that there are no empty entries in the
    # boundaries list and that map correctly yields an inclusion from SRC to TRG
    TRG_final:=List(TRG*1,x->[]);
    map_final:=List(map*1,x->[]);
    for i in [1..Length(TRG)] do
        for j in [1..Length(TRG[i])] do
            if IsBound(TRG[i][j]) then
                if TRG[i][j]<>[0] then
                    if i=1 then
                        Add(TRG_final[1],[1,0]);
                    else
                        cell:=TRG[i][j]*1;
                        for k in [2..Length(cell)] do
                            count:=0;
                            for l in [1..cell[k]] do
                                if not IsBound(TRG[i-1][l]) then
                                    count:=count+1;
                                elif TRG[i-1][l]=[0] then
                                    count:=count+1;
                                fi;
                            od;
                            cell[k]:=cell[k]-count;
                        od;
                        Add(TRG_final[i],cell);
                    fi;
                fi;
            fi;
        od;
    od;
fi;
    return TRG;#TRG_final;
end;
src:=[
    List([1..4],x->[1,0]),
    [
        [2,1,2], [2,1,3], [2,2,4], [2,3,4] 
    ],
    []
];
trg:=[
    List([1..10],x->[1,0]),
    [
        [2,1,2], [2,1,3], [2,2,8], [2,3,4], [2,3,9], [2,4,5], [2,4,6], [2,5,7],
        [2,6,7], [2,7,8], [2,8,10], [2,9,10]
    ],
    [
        [7,1,2,3,4,6,8,10], [7,4,5,7,9,10,11,12], [4,6,7,8,9]
    ],
    []
];
map:=function(i,j)
    if i=0 then
        return j+3;
    elif i=1 then
        return j+5;
    fi;
end;
framed_square:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(src),
        target:=RegularCWComplex(trg),
        mapping:=map
    )
);
src2:=[
    [ [1,0], [1,0] ],
    [ [2,1,2], [2,1,2] ],
    [ ] 
];
trg2:=[
    List([1..4],x->[1,0]),
    [ [2,1,2], [2,1,4], [2,1,4], [2,2,3], [2,2,3], [2,3,4] ],
    [ [4,1,2,4,6], [4,1,3,5,6], [2,4,5], [2,2,3], [2,2,3] ],
    [ [4,1,3,4,5], [4,2,3,4,5] ],
    []
];
map2:=function(i,j)
    if i=0 then
        return j+1;
    else
        return j+3;
    fi;
end;
3ball:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(src2),
        target:=RegularCWComplex(trg2),
        mapping:=map2
    )
);
