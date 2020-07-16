################################################################################
############ Input: an inclusion f: Z -> Y of regular CW-complexes #############
################################################################################
########### Output: let N(Z) denote an open tubular neighbourhood of ###########
################### Z in Y. this will output an inclusion B -> C where C #######
################### is homeomorphic to Y \ N(Z) and B is homeomorphic to the ###
################### boundary of N(Z) ###########################################
################################################################################
RegularCWComplexComplement:=function(f)
    local
        Y, B, IsInternal, path_comp, i,
        j, clsr, int, IsSubpathComponent,
        ext_cell_2_f_notation, f_notation_2_ext_cell,
        ext_cells, k, ext_cell_bnd, e_n_bar,
        ext_cell_cbnd;

    Y:=RegularCWMapToCWSubcomplex(f);
    B:=List([1..Dimension(Y[1])+1],x->[]);

    IsInternal:=function(n,k) # is the kth n-cell of Y in Y\Z?
        if n+1>Length(Y[2]) then
            return true;
        elif k in Y[2][n+1] then
            return false;
        fi;
        return true;
    end;

    path_comp:=List([1..Length(Y[1]!.boundaries)-1],x->[]);
    # list of all of the path components of the intersection between
    # the closure of each internal cell and the subcomplex Z < Y
    # note: entries may be an empty list if they correspond to an
    #       empty intersection or they may not be assigned a
    #       value at all if the associated cell lies in Z
    for i in [1..Length(Y[1]!.boundaries)-1] do
        for j in [1..Length(Y[1]!.boundaries[i])] do
            if IsInternal(i-1,j) then
                Add(B[i],Y[1]!.boundaries[i][j]);
                # output will contain all internal cells
                clsr:=ClosureCWCell(Y[1],i-1,j);
                int:=IntersectionCWSubcomplex(clsr,Y);
                
                # TEST!
                if int[2]<>List(int[2],x->[]) and int[2][2]<>[] then
                    Print(
                        "\n",
                        int,
                        "\n",
                        Homology(Source(CWSubcomplexToRegularCWMap(int)),1),
                        "\n",
                        Homology(Source(CWSubcomplexToRegularCWMap(int)),2),
                        "\n"
                    );
                fi;
                
                path_comp[i][j]:=PathComponentsCWSubcomplex(int);
            else
                Add(B[i],"*"); # temporary entry to keep correct indexing
            fi;
        od;
    od;

    for i in [1..Length(B)] do
        for j in [1..Length(B[i])] do
            if i>1 and B[i][j]<>"*" then
                B[i][j]:=Filtered(
                        B[i][j]{[2..B[i][j][1]+1]},
                        x->"*"<>B[i-1][x]
                    );
                Add(B[i][j],Length(B[i][j]),1);
            fi;
        od;
    od;
    # at this point, B corresponds to the cell complex Y\Z

    IsSubpathComponent:=function(super,sub)
        local
            i;

        if Length(sub[2])<Length(super[2]) then
            for i in [Length(sub[2])+1..Length(super[2])] do
                Add(sub[2],[]);
            od;
        fi;
        for i in [1..Length(super[2])] do
            if not IsSubset(super[2][i],sub[2][i]) then
                return false;
            fi;
        od;

        return true;
    end;

    ext_cell_2_f_notation:=NewDictionary([],true);
    f_notation_2_ext_cell:=NewDictionary([],true);
    # takes a pair [n,k] corresponding to the kth n-cell (external) of B and
    # returns its associated "f notation" i.e. a triple [n+1,k',A] corresponding
    # to the k'th (n+1)-cell (internal) whose closure intersected with Z in the
    # path component A


    ext_cells:=List([1..Length(B)],x->[]);

    for i in [2..Length(path_comp)] do
        for j in [1..Length(path_comp[i])] do
            if IsBound(path_comp[i][j]) then
                if path_comp[i][j]<>[] then
                    # we add as many external (i-1)-cells to B as
                    # there are path components in path_comp[i][j]
                    for k in [1..Length(path_comp[i][j])] do
                        if i=2 then
                            ext_cell_bnd:=[0];
                        else
                            e_n_bar:=ClosureCWCell(Y[1],i-1,j)[2];
                            ext_cell_bnd:=List(
                                ext_cells[i-2],
                                x->LookupDictionary(
                                    ext_cell_2_f_notation,
                                    [i-3,x]
                                )
                            );
                            ext_cell_bnd:=Filtered(
                                ext_cell_bnd,
                                x->
                                    x[2] in e_n_bar[i-1]
                                    and
                                    IsSubpathComponent(path_comp[i][j][k],x[3])
                            );
                            ext_cell_bnd:=List(
                                ext_cell_bnd,
                                x->LookupDictionary(
                                    f_notation_2_ext_cell,
                                    x
                                )[2]
                            );
                        fi;
                        ext_cell_cbnd:=[j];

                        AddCell(
                            B,
                            i-2,
                            ext_cell_bnd,
                            ext_cell_cbnd
                        );
                        Add(ext_cells[i-1],Length(B[i-1]));
                        AddDictionary(
                            ext_cell_2_f_notation,
                            [i-2,Length(B[i-1])],
                            [i-1,j,path_comp[i][j][k]]
                        );
                        f_notation_2_ext_cell!.entries:=List(
                            ext_cell_2_f_notation!.entries,
                            Reversed
                        );
                    od;
                fi;
            fi;
        od;
    od;

    # a final reindexing of B and removal of "*"
    # entries that once corresponded to cells of Z
    for i in [2..Length(B)] do
        for j in [1..Length(B[i])] do
            for k in [2..Length(B[i][j])] do
                B[i][j][k]:=
                B[i][j][k]-
                Length(
                    Filtered(
                        B[i-1]{[1..B[i][j][k]]},
                        x->x="*"
                    )
                );
            od;
        od;
    od;
    B:=List(B,x->Filtered(x,y->y<>"*"));
    Add(B,[]);

    return B;#RegularCWComplex(B);
end
