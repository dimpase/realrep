\\ Inputs: (to be added from GAP...)
\\ cyclotomic_degree = ??? ; desired_norm = ???

\\ Create real and complex cyclotomic fields:
\\ reals = bnfinit(galoissubcyclo(cyclotomic_degree,-1,0,y));
\\ ^^^^^ Supplied from GAP, because we may need it to define the norm.
complexes = rnfinit(reals, x^2 + 1);

\\ Make an element of a given norm:
norms = rnfisnorminit(reals.pol, complexes.pol, 1);
norm_root = rnfisnorm(norms, desired_norm);
\\ return:
cyclotomic_x = liftall(norm_root[1]);
print(cyclotomic_x);