################################################################################
############ Input: a CW-subcomplex [X,S] ######################################
################################################################################
########### Output: a list of CW-subcomplexes [[X,S1],...,[X,Sn]] ##############
################### arising as the path components of [X,S] ####################
################################################################################
PathComponentsCWSubcomplex:=function(XS)
    local   
        ccs, i, j, k, cell, int, l;

    ccs:=List(
        XS[2][1]*1,
        x->Concatenation(
            [[x]],
            List([2..Length(XS[2])],y->[])
        )
    );
    for i in [1..Length(ccs)] do
        for j in [2..Length(ccs[i])] do
            for k in [1..Length(XS[2][j])] do
                cell:=XS[1]!.boundaries[j][XS[2][j][k]]*1;
                cell:=cell{[2..cell[1]+1]};
                int:=Intersection(cell,ccs[i][j-1]);
                if int<>[] then
                    for l in [1..Length(cell)] do
                        Add(ccs[i][j-1],cell[l]);
                    od;
                    Add(ccs[i][j],XS[2][j][k]);
                fi;
            od;
        od;
    od;
    for i in [1..Length(ccs)] do
        for j in Filtered([1..Length(ccs)],x->x<>i) do
            if ccs[i]<>[] and ccs[j]<>[] then
                for k in [1..Length(ccs[i])] do
                    if ccs[j]<>[] then
                        if Intersection(ccs[i][k],ccs[j][k])<>[] then
                            for l in [1..Length(ccs[i])] do
                                ccs[i][l]:=Set(Concatenation(ccs[i][l],ccs[j][l]));
                            od;
                            ccs[j]:=[];
                        fi;
                    fi;
                od;
            fi;
        od;
    od;
    ccs:=Filtered(List(ccs,x->List(x,y->Set(y))),z->z<>[]);

    return List(ccs,x->[ShallowCopy(XS[1]),x]);
end;
