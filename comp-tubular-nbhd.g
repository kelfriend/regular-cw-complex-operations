SubspaceComplement:=function(inc)
    local
        src, trg, dim_s, dim_t, SRC, TRG,
        map, i, j, rep, cell, cobnd, k, x,
        cocell, pos, face, bndbnd, bndbnd_src, l,
        base, int, len, bnd, cc, comp, TRG_final,
        map_final, Filter, count, mapping, ext_faces,
        nest;

    if not IsHapRegularCWMap(inc) then
        Error("the input must be an inclusion of regular CW-complexes");
    fi;

    src:=ShallowCopy(inc!.source);
    trg:=ShallowCopy(inc!.target);

    dim_s:=EvaluateProperty(src,"dimension");
    dim_t:=EvaluateProperty(trg,"dimension");

    SRC:=List([0..dim_t],x->[]);
    TRG:=List(trg!.boundaries*1,x->Concatenation(x,[[0]]));

    map:=List([1..dim_t+1],x->[]);

    # step 1: delete all cells of src from trg
    for i in [0..dim_s] do
        for j in [1..Length(src!.boundaries[i+1])] do
            Unbind(TRG[i+1][inc!.mapping(i,j)]);
        od;
    od;
    ############################################################################

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
    ############################################################################

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
                                            Add(SRC[1],[1,0]);
                                            Add(TRG[1],[1,0]);
                                            Add(map[1],Length(TRG[1]));

                                            Add(int,Length(TRG[1]));
                                            Add(int,2,1);
                                            TRG[2][k]:=int*1;

                                            int:=Concatenation(
                                                [2],
                                                List(int{[2,3]},x->Position(map[1],x))
                                            );
                                            SRC[2][Position(map[2],k)]:=int*1;
                                        elif len=0 then
                                            Add(SRC[1],[1,0]);
                                            Add(TRG[1],[1,0]);
                                            Add(map[1],Length(TRG[1]));

                                            Add(SRC[1],[1,0]);
                                            Add(TRG[1],[1,0]);
                                            Add(map[1],Length(TRG[1]));

                                            int:=[
                                                2,
                                                Length(TRG[1])-1,
                                                Length(TRG[1])
                                            ];
                                            TRG[2][k]:=int*1;

                                            int:=Concatenation(
                                                [2],
                                                List(int{[2,3]},x->Position(map[1],x))
                                            );
                                            SRC[2][Position(map[2],k)]:=int*1;
                                        else
                                            Add(int,2,1);
                                            TRG[2][k]:=int*1;

                                            int:=Concatenation(
                                                [2],
                                                List(int{[2,3]},x->Position(map[1],x))
                                            );
                                            SRC[2][Position(map[2],k)]:=int*1;
                                        fi;
                                    else
                                        Add(int,Length(int),1);
                                        TRG[i-1][k]:=int*1;

                                        int:=Concatenation(
                                            [int[1]],
                                            List(int{[2..Length(int)]},x->Position(map[i-2],x))
                                        );
                                        SRC[i-1][Position(map[i-1],k)]:=int*1;
                                    fi;
                                fi;
                            od; 
                        fi;
                    od;
                fi;
            fi;
        od;
    od;
    ############################################################################

    # step 4: patch all `gaps of type 2', i.e., add cells where necessary in
    # order to ensure that the boundary of the boundary is null in the chain
    # complex
    for i in [3..dim_t+1] do
        for j in [1..Length(TRG[i])*1] do
            if IsBound(TRG[i][j]) then
                if TRG[i][j]<>[] then
                    bndbnd:=TRG[i][j]*1;
                    bndbnd:=bndbnd{[2..bndbnd[1]+1]};
                    bndbnd:=List(bndbnd*1,x->TRG[i-1][x]{[2..TRG[i-1][x][1]+1]});
                    bnd:=List(bndbnd*1,Set);
                    bndbnd:=Concatenation(bndbnd);
                    bndbnd:=Filtered(bndbnd,x->(Length(Positions(bndbnd,x)) mod 2)<>0);
                    bndbnd:=Set(bndbnd);
                    if bndbnd<>[] then # form a (Length(bndbnd))-gon from the edges of TRG[i][j]
                        if i=3 and Length(bndbnd)>2 then
                            cc:=[]; # separate the 0-cells into connected components
                            for k in [1..Length(bnd)] do
                                comp:=[];
                                Add(comp,bnd[k][1]);
                                Add(comp,bnd[k][2]);
                                for l in [k+1..Length(bnd)] do
                                    if comp[Length(comp)] in bnd[l] then
                                        Add(comp,bnd[l][1]);
                                        Add(comp,bnd[l][2]);
                                    fi;
                                od;
                                Add(cc,comp);
                            od;
                            for k in [1..Length(cc)] do
                                if IsBound(cc[k]) then
                                    for l in Difference([1..Length(cc)],[k]) do
                                        if IsBound(cc[l]) then
                                            if Intersection(cc[k],cc[l])<>[] then
                                                cc[k]:=Union(cc[k],cc[l]);
                                                Unbind(cc[l]);
                                            fi;
                                        fi;
                                    od;
                                fi;
                            od;
                            cc:=Set(cc); # now join the ccs together to patch the 2-cell
                            cc:=List(cc,x->Filtered(x,y->y in bndbnd));
                            for k in [1..Length(cc)] do
                                cell:=[
                                    2,
                                    cc[k][2],
                                    cc[(k mod Length(cc))+1][1]
                                ];
                                Add(TRG[2],cell*1);
                                Add(TRG[3][j],Length(TRG[2]));
                                TRG[3][j][1]:=TRG[3][j][1]+1;

                                Add(map[2],Length(TRG[2]));
                                cell:=Concatenation(
                                    [cell[1]],
                                    List(cell{[2..Length(cell)]},x->Position(map[1],x))
                                );
                                Add(SRC[2],cell*1);
                            od;
                        else
                            Add(bndbnd,Length(bndbnd),1);
                            Add(TRG[i-1],bndbnd*1);
                            Add(TRG[i][j],Length(TRG[i-1]));
                            TRG[i][j][1]:=TRG[i][j][1]+1;

                            bndbnd_src:=Concatenation(
                                [bndbnd[1]*1],
                                List(1*bndbnd{[2..Length(bndbnd)]},x->Position(map[i-2],x))
                            );
                            if fail in bndbnd_src then # this shouldn't happen . . .
                                for k in [1..Length(Positions(bndbnd_src,fail))] do
                                    cell:=TRG[i-2][
                                        bndbnd[
                                            Positions(bndbnd_src,fail)[k]
                                        ]
                                    ];
                                    Add(
                                        SRC[i-2],
                                        Concatenation(
                                            [cell[1]],
                                            List(cell{[2..Length(cell)]},x->Position(map[i-3],x))
                                        )
                                    );
                                    Add(map[i-2],bndbnd[Positions(bndbnd_src,fail)[k]]);
                                od;
                            fi;
                            bndbnd_src:=Concatenation(
                                [bndbnd[1]*1],
                                List((bndbnd*1){[2..Length(bndbnd)]},x->Position(map[i-2],x))
                            );
                            Add(map[i-1],Length(TRG[i-1]));
                            Add(SRC[i-1],bndbnd_src*1);
                        fi;
                    fi;
                fi;
            fi;
        od;
    od;
    ############################################################################

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
                        if Position(map[i],j)<>fail then
                            Add(map_final[1],Length(TRG_final[1]));
                        fi;
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
                        
                        if j in map[i] then
                            Add(map_final[i],Length(TRG_final[i]));
                        fi;
                    fi;
                fi;
            fi;
        od;
    od;
    ############################################################################

    # step 5.5: to ensure that the output CW-complex is regular, we must check
    # that each cell added to TRG/TRG_final, g(e^n), is homeomorphic to an open
    # (n-1)-ball
    Filter:=function()
        local
            bool1, map_list, i, j, k,
            bool2, e_n, e_n_bar;

        bool1:=EvaluateProperty(inc,"subspace complement compatibility");

        if bool1<>fail then
            return bool1;
        else
            if not IsBound(inc!.properties) then
                inc!.properties:=[];
            fi;
            map_list:=List([1..Length(trg!.boundaries)],x->[]); # evaluating inc!.mapping
            for i in [1..Length(src!.boundaries)] do # for all possible inputs to save computation
                for j in [1..Length(src!.boundaries[i])] do
                    Add(map_list[i],inc!.mapping(i-1,j));
                od;
            od;
            for i in [3..Length(trg!.boundaries)-1] do
                for j in [1..Length(trg!.boundaries[i])] do
                    if not j in map_list[i] then # we have a cell of X \ Y
                        e_n:=trg!.boundaries[i][j]*1;
                        e_n:=e_n{[2..e_n[1]+1]};
                        e_n_bar:=List([1..i-1],x->[]); # the closure of e_n
                        
                        for k in Reversed([1..i-1]) do
                            for l in [2..Length(e_n)] do
                                Add(e_n_bar[k],e_n[l]);

                                if k>1 then
                                    e_n:=List(
                                        e_n,
                                        x->trg!.boundaries[k][x]{
                                            [2..trg!.boundaries[k][x][1]+1]
                                        }
                                    );
                                    e_n:=Concatenation(e_n);
                                fi;
                            od;
                        od;
                        
                        bool2:=false;
                        for k in [1..i-1] do
                            if Intersection(map_list[k],e_n_bar[k])<>[] then
                                bool2:=true;
                            fi;
                        od;

                        if bool2=false then
                            Add(
                                inc!.properties,
                                ["subspace complement compatibility",false]
                            );
                            return false;
                        fi;
                    fi;
                od;
            od;
            if bool2=true then
                Add(
                    inc!.properties,
                    ["subspace complement compatibility",true]
                );
                return true;
            fi;
        fi;
        
    end;

    #if Filter()=false then
    #    Error("the given complexes are not compatible with this algorithm");
    #fi;
    ############################################################################

    # step 6: adjust SRC so that it is the boundary of TRG with a tubular
    # neighbourhood removed and not just the boundary of the tubular
    # neighbourhood itself
    ext_faces:=[];
    for i in [1..Length(TRG_final[Length(TRG_final)-1])] do
        for j in [2..Length(TRG_final[Length(TRG_final)-1][i])] do
            Add(ext_faces,TRG_final[Length(TRG_final)-1][i][j]);
        od;
    od;
    ext_faces:=Filtered(ext_faces,x->Length(Positions(ext_faces,x))=1);
    ext_faces:=Filtered(ext_faces,x->not x in map_final[Length(map_final)-1]);

    for i in [1..Length(ext_faces)] do
        nest:=List([1..Length(TRG_final)-2],x->[]);
        Add(
            nest[Length(nest)],
            [
                ext_faces[i],
                TRG_final[Length(nest)][ext_faces[i]*1]
            ]
        );
        for j in Reversed([2..Length(nest)]) do
            for k in [1..Length(nest[j])] do
                for l in [2..Length(nest[j][k][2])] do
                    Add(
                        nest[j-1],
                        [
                            nest[j][k][2][l],
                            TRG_final[j-1][nest[j][k][2][l]*1]
                        ]
                    );
                od;
            od;
        od;
        for j in [1..Length(nest)] do
            for k in [1..Length(nest[j])] do
                if not nest[j][k][1] in map_final[j] then
                    Add(map_final[j],nest[j][k][1]);
                    if j>1 then
                        Add(
                            SRC[j],
                            Concatenation(
                                [nest[j][k][2][1]],
                                List(
                                    nest[j][k][2]{[2..nest[j][k][2][1]+1]},
                                    x->Position(map_final[j-1],x)
                                )
                            )
                        );
                    else
                        Add(SRC[1],[1,0]);
                    fi;
                fi;
            od;
        od;
    od;
    ############################################################################

    mapping:={i,j}->map_final[i+1][j];
    
    #return [SRC_final,TRG_final,map_final];
    return Objectify(
        HapRegularCWMap,
        rec(
            source:=RegularCWComplex(SRC),
            target:=RegularCWComplex(TRG_final),
            mapping:=mapping
        )
    );

end;

# this is the first half of the arcknotcomp1s.g code for the example of the unknot
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

# this is the altered decomposition of the 3ball (see 3ball.pdf + email)
src2:=[
    [ [1,0], [1,0] ],
    [ [2,1,2], [2,1,2] ],
    [ ] 
];
trg2:=[
    List([1..4],x->[1,0]),
    [ [2,1,2], [2,1,4], [2,1,4], [2,2,3], [2,2,3], [2,3,4], [2,1,4], [2,1,4], [2,2,3] ],
    [ [4,1,2,4,6], [4,1,3,5,6], [2,4,9], [2,5,9], [2,2,7], [2,3,7], [2,2,8], [2,3,8], [4,1,7,6,9], [4,1,8,6,9] ],
    [ [4,1,9,5,3], [4,2,9,6,4], [4,1,10,7,3], [4,2,10,8,4] ],
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

# this is what 3ball should be if the algorithm works (verified)
src3:=
[ [ [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [1,0], [1,0] ], 
  [ [ 2, 1, 3 ], [ 2, 2, 4 ], [ 2, 1, 3 ], [ 2, 2, 4 ], [ 2, 3, 2 ], 
      [ 2, 4, 1 ], [ 2, 3, 2 ], [ 2, 4, 1 ], [2,5,6], [2,5,6], [2,5,6], [2,5,6] ], 
  [ [ 4, 1, 2, 5, 6 ], [ 4, 3, 4, 5, 6 ], [ 4, 1, 2, 7, 8 ], 
      [ 4, 3, 4, 7, 8 ], [2,9,11],[2,11,10],[2,9,12],[2,10,12] ], [  ] ];
trg3:=
[ [ [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ] ], 
  [ [ 2, 1, 3 ], [ 2, 1, 2 ], [ 2, 1, 2 ], [ 2, 5, 2 ], [ 2, 1, 2 ], 
      [ 2, 1, 2 ], [ 2, 4, 6 ], [ 2, 3, 5 ], [ 2, 4, 6 ], [ 2, 3, 5 ], 
      [ 2, 4, 6 ], [ 2, 5, 4 ], [ 2, 6, 3 ], [ 2, 5, 4 ], [ 2, 6, 3 ] ], 
  [ [ 4, 1, 2, 8, 4 ], [ 4, 1, 3, 10, 4 ], [ 2, 9, 7 ], [ 2, 11, 7 ], 
      [ 2, 2, 5 ], [ 2, 3, 5 ], [ 2, 2, 6 ], [ 2, 3, 6 ], 
      [ 6, 1, 5, 4, 7, 12, 13 ], [ 6, 1, 6, 4, 7, 14, 15 ], 
      [ 4, 8, 9, 12, 13 ], [ 4, 10, 11, 12, 13 ], [ 4, 8, 9, 14, 15 ], 
      [ 4, 10, 11, 14, 15 ] ], 
  [ [ 5, 1, 9, 5, 3, 11 ], [ 5, 2, 9, 6, 4, 12 ], [ 5, 1, 10, 7, 3, 13 ], 
      [ 5, 2, 10, 8, 4, 14 ] ], [  ] ];
mp3:=
[ [ 3, 4, 5, 6, 1, 2 ], [ 8, 9, 10, 11, 12, 13, 14, 15, 2,3,5,6 ], [ 11, 12, 13, 14 , 5,6,7,8], ];
map3:={i,j}->mp3[i+1][j];
test:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(src3),
        target:=RegularCWComplex(trg3),
        mapping:=map3
    )
);

# this is a standard cube with a subcomplex consisting of half of the edges (see cube2.pdf)
src4:=
[
    [
        List([1..7],x->[1,0]),
        [
            [2,1,2],[2,2,3],[2,3,4],[2,4,5],[2,5,6],[2,6,7]
        ],
        [

        ]
    ]
];
trg4:=
[
    List([1..8],x->[1,0]),
    [
        [2,1,2],[2,1,4],[2,1,5],[2,2,3],[2,2,6],[2,3,4],[2,3,7],[2,4,8],[2,5,6],[2,5,8],[2,6,7],[2,7,8]
    ],
    [
        [4,1,2,4,6],[4,2,3,8,10],[4,1,3,5,9],[4,4,5,7,11],[4,6,5,8,12],[4,9,10,11,12]
    ],
    [
        [6,1,2,3,4,5,6]
    ],
    [

    ]
];
mp4:=[[6,5,1,2,3,7,8],[9,3,1,4,7,12]];
map4:={i,j}->mp4[i+1][j];
cube:=Objectify(
    HapRegularCWMap,
    rec(
        source:=RegularCWComplex(src4),
        target:=RegularCWComplex(trg4),
        mapping:=map4
    )
);