
\\ The following function computes the Gauss-Manin connection for a family of hyperelliptic curves
\\ y^2 = P(x,t)

GaussManin(P)={

    Px = deriv(P,x);
    Pt = deriv(P,t);

    [A,B,C] = polresultantext(P,Px);

    A = A/C;
    B = B/C;

    At = deriv(A,t);
    Bx = deriv(B,x);
    Bt = deriv(B,t);


    \\ AP+BPx=1

    \\ omega = Aydx + 2Bdy

    domega = -1/2*A*Pt - At*P + Bx*Pt-Bt*Px;

    dxomega = x*domega + B*Pt;

    deg_domega = if(domega == 0,0,poldegree(domega,x));
    deg_dxomega = poldegree(dxomega,x);

    deg_max = max(deg_domega,deg_dxomega);

    vec_domega = if(domega == 0,[0],Vec(domega));
    vec_dxomega = Vec(dxomega);

    vec_P = Vec(P);
    vec_Px = Vec(Px);

    vec_coefx = Vec(0,deg_max+1);
    vec_coefx[deg_max] = 1;
    vec_coef0 = Vec(0,deg_max+1);
    vec_coef0[deg_max+1] = 1;

    for(n=1,deg_max-1,

        M_P = matrix(n,n,i,j,if(i==j,vec_Px[1],if(i==j+1,vec_Px[2],if(i==j+2,vec_Px[3],0))));
        M_Px = matrix(n,n,i,j,if(i==j,(n-i)*vec_P[1],if(i==j+1,(n-i+1)*vec_P[2],if(i==j+2,(n-i+2)*vec_P[3],if(i==j+3,(n-i+3)*vec_P[4],0)))));
        
        M = 1/2*M_P+M_Px;

        unit_vec_xn = Vec(1,n)~;

        F = Pol(M^(-1)*unit_vec_xn);

        xn_as_exact_form = 1/2*F*deriv(P,x)+deriv(F,x)*P;
        \\ xn_as_exact_form =  X^(n+1)+coeff_x*X+coeff_0

        coeff_x = Vec(xn_as_exact_form)[n+2-1];
        vec_coefx[deg_max-n] = -1*coeff_x;

        coeff_0 = Vec(xn_as_exact_form)[n+2];
        vec_coef0[deg_max-n] = -1*coeff_0;

    );

    xomega_coeff_domega= vec_coefx[(deg_max-deg_domega+1)..deg_max+1]*vec_domega~;
    omega_coeff_domega= vec_coef0[(deg_max-deg_domega+1)..deg_max+1]*vec_domega~;

    xomega_coeff_dxomega= vec_coefx*vec_dxomega~;
    omega_coeff_dxomega= vec_coef0*vec_dxomega~;

    GMmatrix = [omega_coeff_domega,xomega_coeff_domega;omega_coeff_dxomega,xomega_coeff_dxomega];
    
    monodromy = Monodromy(GMmatrix);
    
    return(GMmatrix);

}

Monodromy(GMmatrix)={

    residue_matrix = matrix(#GMmatrix,#GMmatrix,i,j,subst(GMmatrix[i,j]*t,t,0));
    monodromy = matrix(#GMmatrix);

    for(i=0,100,monodromy+= (2*Pi*I*residue_matrix)^i/factorial(i));

    return(monodromy);
}


Period_Matrix(e1,e2,e3)={

   \\ given y^2 = (x-e1)(x-e2)(x-e3) we want to evaluate the periods e1-e2 and e1-e3

    lambda = (e3-e1)/(e2-e1);

    \\ period matrix Legendre family w/ 0-1 0-lambda periods
    Period_M = Pi*[1/sqrt(lambda)*hypergeom([1/2,1/2],1,1/lambda),hypergeom([1/2,1/2],1,lambda); 
                    1/(2*sqrt(lambda))*hypergeom([3/2,1/2],2,1/lambda),lambda/2*hypergeom([3/2,1/2],2,lambda)];

    \\base change sending e1 - 0, e_2 - 1, e3 - lambda
    Period_M = [1/sqrt(e2-e1),0;
                e1/sqrt(e2-e1),sqrt(e2-e1)]*Period_M*[0,1;1,0];

    
    return(Period_M);

}


exponential_matrix(M)={

    exp_m = 0;
    for(i=0,100,exp_m+= M^i/factorial(i));

    return(exp_m);

}

GaussManinHyp(P)={

    deg_P = poldegree(P,x);
    Px = deriv(P,x);
    Pt = deriv(P,t);

    [A,B,C] = polresultantext(P,Px);

    A = A/C;
    B = B/C;

    At = deriv(A,t);
    Bx = deriv(B,x);
    Bt = deriv(B,t);


    \\ AP+BPx=1

    \\ omega = Aydx + 2Bdy

    domega = -1/2*A*Pt - At*P + Bx*Pt-Bt*Px;

    differentials = matrix(1,deg_P-1,i,j,domega*x^(j-1)+(j-1)*x^(j-2)*B*Pt)[1,];

    \\dxomega = x*domega + B*Pt;
    degree_differentials = matrix(1,deg_P-1,i,j,poldegree(differentials[j],x))[1,];
    
    \\deg_domega = if(domega == 0,0,poldegree(domega,x));
    \\deg_dxomega = poldegree(dxomega,x);

    deg_max = vecmax(degree_differentials);

    vec_differentials = matrix(1,deg_P-1,i,j,Vec(differentials[j]))[1,];
    
    \\vec_domega = if(domega == 0,[0],Vec(domega));
    \\vec_dxomega = Vec(dxomega);

    vec_P = Vec(P);
    vec_Px = Vec(Px);

    vec_coefx = Vec(0,deg_max+1);
    vec_coefx[deg_max] = 1;
    vec_coef0 = Vec(0,deg_max+1);
    vec_coef0[deg_max+1] = 1;

    coeff_xn = matrix(deg_max-1+deg_P-1,deg_P-1,i,j,if(i>(deg_max-1) && (i-deg_max+1)==j,1));


    for(n=1,deg_max-1,

        \\M_P = matrix(n,n,i,j,if(i==j,vec_Px[1],if(i==j+1,vec_Px[2],if(i==j+2,vec_Px[3],0))));
        M_P = matrix(n,n,i,j,if(i-j>=0 && i-j<#vec_Px,vec_Px[i-j+1]));
        \\M_Px = matrix(n,n,i,j,if(i==j,(n-i)*vec_P[1],if(i==j+1,(n-i+1)*vec_P[2],if(i==j+2,(n-i+2)*vec_P[3],if(i==j+3,(n-i+3)*vec_P[4],0)))));
        M_Px = matrix(n,n,i,j,if(i-j>=0 && i-j<#vec_P,(n-j)*vec_P[i-j+1]));

        M = 1/2*M_P+M_Px;

        unit_vec_xn = Vec(1,n)~;

        F = Pol(M^(-1)*unit_vec_xn);

        xn_as_exact_form = 1/2*F*deriv(P,x)+deriv(F,x)*P;
        \\ xn_as_exact_form =  X^(n+1)+coeff_x*X+coeff_0
        
        coeff_xn[deg_max-n,] = (-1)*Vec(xn_as_exact_form)[n+1..n+deg_P-1];

        \\coeff_x = Vec(xn_as_exact_form)[n+2-1];
        \\vec_coefx[deg_max-n] = -1*coeff_x;

        \\coeff_0 = Vec(xn_as_exact_form)[n+2];
        \\vec_coef0[deg_max-n] = -1*coeff_0;

    );

    \\xomega_coeff_domega= vec_coefx[(deg_max-deg_domega+1)..deg_max+1]*vec_domega~;
    \\omega_coeff_domega= vec_coef0[(deg_max-deg_domega+1)..deg_max+1]*vec_domega~;

    \\xomega_coeff_dxomega= vec_coefx*vec_dxomega~;
    \\omega_coeff_dxomega= vec_coef0*vec_dxomega~;

    GMmatrix = matrix(deg_P-1,deg_P-1,i,j,((coeff_xn[,deg_P-j])~)[deg_max-1+deg_P-#vec_differentials[i]..deg_max-1+deg_P-1]*(vec_differentials[i])~);

    \\GMmatrix = [omega_coeff_domega,xomega_coeff_domega;omega_coeff_dxomega,xomega_coeff_dxomega];
    
    \\monodromy = Monodromy(GMmatrix);
    
    return(GMmatrix);

}

Frobenius(P,p)={

    deg_P = poldegree(P,x);
    Px = deriv(P,x);
    Pt = deriv(P,t);

    [A,B,C] = polresultantext(P,Px);

    A = A/C;
    B = B/C;

    At = deriv(A,t);
    Bx = deriv(B,x);
    Bt = deriv(B,t);

    \\ AP+BPx=1

    \\ omega = Aydx + 2Bdy

    \\F_P^*(omega)=x^(p-1)F_p(y)^(-1)dx

    F_pz =y^(-p)*(1+y^(-2*p)*(subst(P,x,x^p)-P^p))^(-1/2); 

    domega = -1/2*A*Pt - At*P + Bx*Pt-Bt*Px;

    Fp_differentials = matrix(1,deg_P-1,i,j,p*x^(p*(j-1)+(p-1))*F_pz)[1,];

    \\dxomega = x*domega + B*Pt;
    \\degree_differentials = matrix(1,deg_P-1,i,j,poldegree(Fp_differentials[j],x))[1,];
    
}