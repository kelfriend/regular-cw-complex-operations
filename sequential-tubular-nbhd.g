SequentialTubularNeighbourhood:=function(inc)
    local
        sub, Y, tub_sub, i, closure,
        j, k, bnd, l, pos, cell, reind,
        max;

    sub:=RegularCWMapToCWSubcomplex(inc);
    Y:=sub[1]!.boundaries*1;
    sub:=sub[2]*1;
    tub_sub:=List([1..Length(sub)+1],x->[]);

    # loop through all cells of Y of dimension dim(Y)
    for i in [1..Length(Last(sub))] do
    # first find the closure of the Last(sub)[i]th (Length(sub)-1)-cell
        closure:=List([1..Length(sub)],x->[]);
        Add(Last(closure),Last(sub)[i]);
        for j in Reversed([2..Length(sub)]) do
            for k in [1..Length(closure[j])] do
                bnd:=[];
                for l in [2..Length(Y[j][closure[j][k]])] do
                    Add(bnd,Y[j][closure[j][k]][l]);
                od;
                closure[j-1]:=Union(closure[j-1],bnd);
            od;
        od;
        # remove everything in the closure from Y
        for j in [1..Length(closure)] do
            for k in [1..Length(closure[j])] do
                Y[j][closure[j][k]]:=[0];
                # now unbind any removed cells from
                # the cells in their former coboundaries
                for l in [1..Length(Y[j+1])] do
                    pos:=Position(
                        Y[j+1][l]{
                            [2..Y[j+1][l][1]+1]
                        },
                        closure[j][k]
                    );
                    if pos<>fail then
                        cell:=Y[j+1][l]*1;
                        cell:=cell{
                            Difference(
                                [2..cell[1]+1],
                                [pos+1]
                            )
                        };
                        Add(cell,Length(cell),1);
                        Y[j+1][l]:=cell*1;
                    fi;
                od;
                if closure[j][k] in tub_sub[j] then
                    Unbind(tub_sub[j][Position(tub_sub[j],closure[j][k])]);
                    tub_sub[j]:=Set(tub_sub[j]);
                fi;
            od;
        od;
        # `repair' the boundaries of all cells
        # effected by the above deletion by adding new
        # cells where needed
        for j in [2..Length(Y)] do
            for k in [1..Length(Y[j])] do
                if j=2 then
                    if Y[j][k][1]=1 then
                        AddCell(Y,0,[0],[k]);
                        Add(tub_sub[j-1],Length(Y[j-1]));
                    fi;
                else
                    bnd:=Y[j][k]{
                        [2..Y[j][k][1]+1]
                    };
                    bnd:=List(
                        bnd,
                        x->Y[j-1][x]{
                            [2..Y[j-1][x][1]+1]
                        }
                    );
                    bnd:=Concatenation(bnd);
                    bnd:=Filtered(bnd,x->(Length(Positions(bnd,x)) mod 2)<>0);
                    if bnd<>[] then
                        AddCell(Y,j-2,bnd,[k]);
                        Add(tub_sub[j-1],Length(Y[j-1]));
                    fi;
                fi;
            od;
        od;
    od;

    # remove all [0] entries and reindex everything
    reind:=List([1..Length(Y)-1],x->Positions(Y[x],[0]));
    reind:=List(reind,x->Filtered(x,y->not y+1 in x));
    for i in [1..Length(reind)] do
        for j in [1..Length(reind[i])] do
            reind[i][j]:=[
                reind[i][j],
                Length(
                    Filtered(
                        Y[i]{[1..reind[i][j]]},
                        x->x=[0]
                    )
                )
            ];
        od;
    od;

    Y:=List(Y,x->Filtered(x,y->y<>[0]));
    for i in [2..Length(Y)-1] do
        if reind[i-1]<>[] then
            for j in [1..Length(Y[i])] do
                for k in [2..Length(Y[i][j])] do
                    max:=Filtered(
                        List(reind[i-1],x->x[1]),
                        y->Y[i][j][k]>y
                    );
                    if max<>[] then
                        max:=Maximum(max);
                        max:=reind[i-1][
                            Position(
                                List(reind[i-1],x->x[1]),
                                max
                            )
                        ][2];
                        Y[i][j][k]:=Y[i][j][k]-max;
                    fi;
                od;
            od;
        fi;
    od;
    # repeat the above for tub_sub 
    # (this is done separately to avoid reindexing multiple times)
    for i in [1..Length(tub_sub)] do
        if reind[i]<>[] then
            for j in [1..Length(tub_sub[i])] do
                max:=Filtered(
                        List(reind[i],x->x[1]),
                        y->tub_sub[i][j]>y
                    );
                if max<>[] then
                    max:=Maximum(max);
                    max:=reind[i][
                        Position(
                            List(reind[i],x->x[1]),
                            max
                        )
                    ][2];
                    tub_sub[i][j]:=tub_sub[i][j]-max;
                fi;
            od;
        fi;
    od;

    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(Y),
            tub_sub
        ]
    );
end;
