# BarycentricallySubdivideCell
################################################################################
############ Input: an inclusion of cell complexes f:Z->Y and two integers #####
################### k>0 & n>=0 #################################################
################################################################################
########### Output: the regular CW-map f':Z'->Y' corresponding to the same #####
################### spaces but where the kth n-cell of Y has been ##############
################### barycentrically subdivided #################################
################################################################################
BarycentricallySubdivideCell:=function(f,n,k)
    local
        inc, Y, clsr, complex,
        Subdivide, i, j;

    inc:=RegularCWMapToCWSubcomplex(ShallowCopy(f));
    Y:=inc[1];

    clsr:=ClosureCWCell(Y,n,k);
    complex:=Y!.boundaries*1;
    clsr:=clsr[2]*1;

    Subdivide:=function(a,b)
    # subdivides the bth a-cell into as many a-cells
    # as there are in its barycentric subdivision
        local
            sub_clsr, i, j, l, ints, bnd_ints,
            bnd, pre_len, m, pre_len_bnd;

        sub_clsr:=ClosureCWCell(RegularCWComplex(complex),a,b)[2];

        Add(complex[1],[1,0]); # the barycentre
        if Length(inc[2])>=a+1 then
            if b in inc[2][a+1] then
                Add(inc[2][1],Length(complex[1]));
            fi;
        fi;

        pre_len_bnd:=List([1..a],x->Length(complex[x]));

        for i in [1..a] do
            pre_len:=Length(complex[i+1])*1;
            for j in [1..Length(sub_clsr[i])] do
                if i=1 then # edge case is dealt with separately
                    if a=1 and j=1 then # we overwrite the bth a-cell with
                    # the first a-cell of its barycentric subdivision
                        complex[2][b]:=[
                            2,
                            sub_clsr[1][1],
                            Length(complex[1])
                        ];
                    else
                        Add(
                            complex[2],
                            [
                                2,
                                sub_clsr[1][j],
                                Length(complex[1])
                            ]
                        );
                        if Length(inc[2])>=a+1 then
                            if b in inc[2][a+1] then
                                Add(inc[2][2],Length(complex[2]));
                            fi;
                        fi;
                    fi;
                else
                    bnd:=[];
                    ints:=List([pre_len_bnd[i]+1..Length(complex[i])]);
                    bnd_ints:=List(
                        ints,
                        x->complex[i][x]{[2..complex[i][x][1]+1]}
                    );
                    bnd_ints:=List(
                        bnd_ints,
                        x->Intersection(
                            x,
                            complex[i][sub_clsr[i][j]]{
                                [2..complex[i][sub_clsr[i][j]][1]+1]
                                }
                        )<>[]
                    );
                    bnd:=Concatenation(
                        [sub_clsr[i][j]],
                        Filtered(ints,x->bnd_ints[Position(ints,x)]=true)
                    );
                    Add(bnd,Length(bnd),1);
                    if i=a and j=1 then # overwrite as before
                        complex[a+1][b]:=bnd;
                    else
                        Add(complex[i+1],bnd);

                        if Length(inc[2])>=a+1 then
                            if b in inc[2][a+1] then
                                Add(inc[2][i+1],Length(complex[i+1]));
                            fi;
                        fi;
                    fi;
                fi;

                if i=a and a<Dimension(Y) then # rewrite the coboundary of
                # the replaced a-cell to contain all the a-cells in its
                # barycentric subdivision
                    for l in [2..Length(Y!.coboundaries[a+1][b])] do
                        for m in [pre_len+1..Length(complex[a+1])] do
                            Add(
                                complex[a+2][Y!.coboundaries[a+1][b][l]],
                                m
                            );
                        od;
                        complex[a+2][Y!.coboundaries[a+1][b][l]][1]:=Length(
                            complex[a+2][Y!.coboundaries[a+1][b][l]]
                        )-1;
                    od;
                fi;
            od;
        od;  
    end;

    for i in [2..Length(clsr)] do # inductively apply Subdivide to obtain the
        for j in [1..Length(clsr[i])] do # barycentric subdivision of the
            Subdivide(i-1,clsr[i][j]); # kth n-cell of Y as required
        od;
    od;

    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(complex),
            inc[2]
        ]
    );
end;
