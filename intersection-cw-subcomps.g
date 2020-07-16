################################################################################
############ Input: two CW-subcomplexes [X,S], [X,S'] ##########################
################################################################################
########### Output: the CW-subcomplex [X,S''] corresponding to #################
################### their intersection #########################################
################################################################################
IntersectionCWSubcomplex:=function(XS_1,XS_2)
    local
        max, xs_1, xs_2, i, j;
    
    max:=Maximum(Length(XS_1[2]),Length(XS_2[2]));
    xs_1:=ShallowCopy(XS_1[2]);
    xs_2:=ShallowCopy(XS_2[2]);

    for i in [xs_1,xs_2] do
        if Length(i)<max then
            for j in [Length(i)+1..max] do
                Add(i,[]);
            od;
        fi;
    od;

    return[
        ShallowCopy(XS_1[1]),
        List(
            [1..max],
            x->Intersection(xs_1[x],xs_2[x])
        )
    ];
end;
