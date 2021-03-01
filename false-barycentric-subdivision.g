SubdivideCell:=function(f,n,k)
    local
        sub, Y, closure, plus1,
        i, j, bnd;

    sub:=RegularCWMapToCWSubcomplex(ShallowCopy(f));
    Y:=sub[1]; # the actual CW-complex
    sub:=sub[2]*1; # the indexing of the subcomplex
    closure:=ClosureCWCell(Y,n,k)[2];
    Y:=Y!.boundaries*1;

    Add(Y[1],[1,0]); # the barycentre of the kth n-cell
    plus1:=List([1..Length(closure)],x->[]);
    # this will associate an x-cell in the closure to
    # the resulting (x+1)-cell in the subdivision

    for i in [1..Length(closure)-1] do
        for j in [1..Length(closure[i])] do
            if i=1 then
                Add(Y[2],[2,closure[i][j],Length(Y[1])]);
                Add(plus1[1],Length(Y[2]));
            else
                bnd:=Y[i][closure[i][j]];
                bnd:=bnd{[2..bnd[1]+1]};
                bnd:=List(bnd,x->plus1[i-1][Position(closure[i-1],x)]);
                Add(bnd,closure[i][j],1);
                Add(bnd,Length(bnd),1);
                if i=Length(closure)-1 and j=1 then
                    Y[i+1][closure[i+1][1]]:=bnd;
                    Add(plus1[i],closure[i+1][1]);
                else
                    Add(Y[i+1],bnd);
                    Add(plus1[i],Length(Y[i+1]));
                fi;
            fi;
        od;
    od;
    for i in [1..Length(Y[Length(closure)])] do
        if Last(closure)[1] in
            Y[Length(closure)][i]{[2..Y[Length(closure)][i][1]+1]} then
                Append(Y[Length(closure)][i],plus1[Length(closure)]);
                Unbind(Y[Length(closure)][i][1]);
                Y[Length(closure)][i]:=Set(Y[Length(closure)][i]);
                Add(Y[Length(closure)][i],Length(Y[Length(closure)][i]),1);
        fi;
    od;
    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(Y),
            sub
        ]
    );
end;
