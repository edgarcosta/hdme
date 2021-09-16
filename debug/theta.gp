
/* ---------- Computing with theta characteristics ---------- */

\\ Compute all 10 even Theta characteristics
theta_evenchars() =
{
   my(c1, c2, c3, c4, T);
   T = List([]);
   for(c1=0,1,
      for(c2=0,1,
	 for(c3=0,1,
	    for(c4=0,1,
	       if ([c1,c2]*[c3,c4]~ %2 == 0,
	       listput(T, [[c1, c2], [c3, c4]]));
	    );
	 );
      );
   );
   T;
};

{addhelp(theta_evenchars, "theta_evenchars(): return the list of all
	 10 even theta characteristics.");};

\\ Compute all syzygous triples of theta characteristics
\\ Here T = theta_evenchars()
theta_is_goepel(a1, a2, a3, a4) = {(a1+a2+a3+a4) %2 == 0};

{addhelp(theta_is_goepel, "theta_is_goepel(a1, a2, a3, a4): decide
	 whether [a1, a2, a3, a4] is a Goepel quadruple of theta
	 characteristics.");};

theta_goepel_quadruples(T) =
{
    my(res, i, j, k, l, a1, a2, a3, a4);
    res = List([]);
    for (i=1, 10,
         for (j=i+1, 10,
	      for (k=j+1, 10,
	           for (l=k+1, 10,
		        a1 = T[i];
		        a2 = T[j];
		        a3 = T[k];
			a4 = T[l];
		        if (theta_is_goepel(a1, a2, a3, a4),
			    listput(res, [a1,a2,a3,a4]);
			   );
		       );
		  );
             );
	);
    res;
}

{addhelp(theta_goepel_quadruples, "theta_goepel_quadruples(T): given T
	 = theta_evenchars(), return the list of all Goepel quadruples
	 of theta characteristics.");};

theta_is_syzygous(a1, a2, a3, T) =
{
   my(i, s);
   if (a1 == a2 || a2 == a3 || a3 == a1,
       return(0);
      );
   for (i = 1, 10,
	s = T[i];
	if (s != a1 && s != a2 && s != a3 && theta_is_goepel(a1,a2,a3,s),
	    return(1);
	   );
       );
   return(0);
}

{addhelp(theta_is_syzygous, "theta_is_syzygous(a1, a2, a3, T): given T
	 = theta_evenchars(), decide whether [a1, a2, a3] is a
	 syzygous triple of theta characteristics.");};

theta_syzygous_triples(T) =
{
   my(res, i, j, k, a1, a2, a3);
   res = List([]);
   for (i = 1, 10,
	for (j = i+1, 10,
	     for (k = j+1, 10,
		  a1 = T[i];
		  a2 = T[j];
		  a3 = T[k];
		  if (theta_is_syzygous(a1, a2, a3, T),
		      listput(res, [a1, a2, a3]);
		     );
		 );
	    );
       );
   res;
}

{addhelp(theta_syzygous_triples, "theta_syzygous_triples(T): given T =
	 theta_evenchars(), return the list of all syzygous triples of
	 theta characteristics.");};


/* ---------- Dupont and Streng numerotation for even theta constants ---------- */

\\ Compute the theta index corresponding to the characteristic
theta_label_from_char(a) =
{
   a[2][1] + 2*a[2][2] + 4*a[1][1] + 8*a[1][2];
}

{addhelp(theta_label_from_char, "theta_label_from_char(a): compute the Dupont (or Streng)
	 label (between 0 and 15) corresponding to the theta
	 characteristic a".);};

theta_char_from_label(l) =
{
   my(a, r);
   a = [[0,0],[0,0]];
   [l, r] = divrem(l, 2);
   a[2][1] = r;
   [l, r] = divrem(l, 2);
   a[2][2] = r;
   [l, r] = divrem(l, 2);
   a[1][1] = r;
   [l, r] = divrem(l, 2);
   a[1][2] = r;
   if (a[1] * a[2]~ % 2 != 0, error("Theta constant number ", l, " is not even"););
   a;    
}

theta_vecindex_from_label(l) =
{
   \\{0,1,2,3,4,6,8,9,12,15}+1
   [1,2,3,4,5,0,6,0,7,8,0,0,9,0,0,10][l+1];
}

{addhelp(theta_vecindex_from_label, "theta_vecindex_from_label(l):
	 given a Dupont (or Streng) label l such that the
	 corresponding theta constant is even, return its index in the
	 10theta2 vector. Otherwise, return 0.");};

theta_label_from_vecindex(i) =
{
   [0,1,2,3,4,6,8,9,12,15][i];
};

\\ Returns an error if a is not an even theta characteristic.
get_theta(a, v) =
{
   v[theta_vecindex_from_label(theta_label_from_char(a))];
}

{addhelp(get_theta, "get_theta(a, v): given an even theta
	 characteristic a, and v the vector output by cmh_10theta2_*,
	 return the corresponding theta constant.");};

/* ---------- Formulae for Igusa invariants ---------- */

theta_h10(v) =
{
   vecprod(v);
}

{addhelp(theta_h10, "theta_h10(v): given the vector v containing the
	 squares of the 10 even theta constants at tau, compute
	 h10(tau) in Streng's notation.");};

theta_h4(v) =
{
   my(res, T, i);
   res = 0;
   T = theta_evenchars(); \\there are 10 of them.
   for (i = 1, 10,
	res += get_theta(T[i], v)^4;
       );
   res;
}

{addhelp(theta_h4, "theta_h4(v): given the vector v containing the
	 squares of the 10 even theta constants at tau, compute
	 h4(tau) in Streng's notation.");};

\\Compute the sign in the expression for h6
theta_h6_pm(b, c, d) =
{
   my(b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,d3,d4, e);
   b1 = b[1][1];
   b2 = b[1][2];
   b3 = b[2][1];
   b4 = b[2][2];
   c1 = c[1][1];
   c2 = c[1][2];
   c3 = c[2][1];
   c4 = c[2][2];
   d1 = d[1][1];
   d2 = d[1][2];
   d3 = d[2][1];
   d4 = d[2][2];
   e = b1 + b2 + c1 + c2 + d1 + d2 + b1*c1 + b2*c2 + b4*c2 + b1*c3 - b2*c4 +
      b1*d1 - b3*d1 + c1*d1 + b2*d2 + c2*d2 + c4*d2 + c1*d3 - b2*b3*c1 -
      b2*b4*c2 - b1*b2*c3 - b2*b3*d1 - b3*c1*d1 - b1*c3*d1 - b2*c3*d1
      - b2*b4*d2 - b4*c2*d2 - b1*b2*d3 - b1*c1*d3 - b2*c1*d3;
   (-1)^e;
}

{addhelp(theta_h6_pm, "theta_h6_pm(b, c, d): given a syzygous triple
	 [b,c,d] of even theta characteristics, return the
	 corresponding sign (+-1) in the formula for h6.");};

theta_h6(v) =
{
   my(res, T, triples, i, b, c, d, sgn, term);
   res = 0;
   T = theta_evenchars();
   triples = theta_syzygous_triples(T); \\60 of them.
   for (i = 1, 60,
	[b, c, d] = triples[i];
	sgn = theta_h6_pm(b, c, d);
	term = get_theta(b, v) * get_theta(c, v)
	      * get_theta(d, v);
	res += sgn * term^2;
       );
   res;
}

{addhelp(theta_h6, "theta_h6(v): given the vector v containing the
	 squares of the 10 even theta constants at tau, compute
	 h6(tau) in Streng's notation.");};

theta_h12(v) =
{
    my(res, term, T, S, i, C, j, c);
    T = theta_evenchars();
    S = theta_goepel_quadruples(T);
    res = 0;
    for (i=1, length(S),
         term = 1;
         C = S[i];
	 for (j = 1, 10,
	      c = T[j];
	      if (c != C[1] && c != C[2] && c != C[3] && c != C[4],
		  term *= get_theta(c, v);
		 );
	     );
	 res += term^2;
	);
    res;
}

{addhelp(theta_h12, "theta_h12(v): given the vector v containing the
	 squares of the 10 even theta constants at tau, compute
	 h12(tau) in Streng's notation.");};

theta_h16(v) =
{
    my(res, term1, term2, T, S, i, C, j, c);
    T = theta_evenchars();
    S = theta_goepel_quadruples(T);
    res = 0;
    for (i=1, length(S),
	 term1 = 1;
	 term2 = 0;
	 C = S[i];
	 \\Compute term1 = \prod_{c in T\C} theta[c]^2 as in h12
	 for (j = 1, 10,
	      c = T[j];
	      if (c != C[1] && c != C[2] && c != C[3] && c != C[4],
		  term1 *= get_theta(c, v);
		 );
	     );
	 \\Compute term2 = \sum_{d in C} theta[d]^8
	 for (j = 1, 4,
	      c = C[j];
	      term2 += get_theta(c,v)^4;
	     );     
	 res += term2 * term1^2;
	);
    res;
}

{addhelp(theta_h16, "theta_h16(v): given the vector v containing the
	 squares of the 10 even theta constants at tau, compute
	 h16(tau) in Streng's notation.");};

h4h6h10h12_from_theta(v) =
{
   [theta_h4(v),
    theta_h6(v),
    theta_h10(v),
    theta_h12(v)];
};

{addhelp(h4h6h10h12_from_theta, "h4h6h10h12_from_theta(v): given the
	 vector v containing the squares of the 10 even theta
	 constants at tau, compute [h4(tau), h6(tau), h10(tau),
	 h12(tau)].");};

\\ Consistent with the builtin cmh_I2I4I6I10 function.
igusa_from_theta(v) =
{
   my(h4,h6,h10,h12,j1,j2,j3);
   h4 = theta_h4(v);
   h6 = theta_h6(v);
   h10 = theta_h10(v);
   h12 = theta_h12(v);
   \\ I4 I6'/I10 = h4 h6/h10
   j1 = h4 * h6 / h10;
   \\ I2 I4^2/I10 = h12 h4^2/h10^2
   j2 = h12 * h4^2 / h10^2;
   \\ I4^5 / I10^2 = h4^5 / h10^2
   j3 = h4^5 / h10^2;
   [j1,j2,j3];
};

{addhelp(igusa_from_theta, "igusa_from_theta(v): given the vector v
	 containing the squares of the 10 even theta constants at tau,
	 compute the Igusa--Streng invariants [j1(tau), j2(tau),
	 j3(tau)].");};
