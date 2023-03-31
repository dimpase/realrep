# Adapted from https://github.com/gap-packages/alnuth/blob/e13a7c36159c1475d0d2d63b9dea54fb36873fec/gap/kantin.gi

DeclareGlobalFunction( "PariNorms" );
InstallGlobalFunction(PariNorms, function(cyclotomic_degree, norm_mu)
    local conversion_matrix, inputs, output, paricode, parsed_output, power, real_degree, subfield_coeffs, symmetric_roots, x, y;

    # Test, whether AL_EXECUTABLE is set:
    if AL_EXECUTABLE = fail then
        Error( "AL_EXECUTABLE, the executable for PARI/GP, has to be set" );
    fi;

    # add the prepared code fragments for the calculations in PARI/GP
    # paricode := InputTextFile(Filename(AL_PATH, "norm_equation.gp"));
    paricode := InputTextFile(Filename(DirectoryCurrent(), "norm_equation.gp"));

    # PARI/GP inputs *** Does GAP have something like f-strings? ***:
    Print(norm_mu, "\n\n");
    inputs := ["cyclotomic_degree = "];
    Add(inputs, String(cyclotomic_degree));
    # First we make the real cyclotomic subfield, so that norm is always expressible.
    Add(inputs, ";\nreals = bnfinit(galoissubcyclo(cyclotomic_degree,-1,0,y));\n");
    inputs := Concatenation(inputs, ["desired_norm = "]);
    if norm_mu in Rationals then  # Easy case
        Add(inputs, String(norm_mu));
        Add(inputs, ";\n");
    else  # Need to parse the E(n)-expression into GP-readable form:
        real_degree := Phi(cyclotomic_degree) / 2;  # Always integer.
        symmetric_roots := E(cyclotomic_degree) + ComplexConjugate(E(cyclotomic_degree));
        conversion_matrix := [];
        for power in [0..real_degree] do
            Add(conversion_matrix, CoeffsCyc(symmetric_roots ^ power, cyclotomic_degree));
        od;
        subfield_coeffs := SolutionMat(conversion_matrix, CoeffsCyc(norm_mu, cyclotomic_degree));
        # Print(subfield_coeffs);
        for power in [0..real_degree] do
            Add(inputs, String(subfield_coeffs[power + 1]));
            Add(inputs, "*y^");
            Add(inputs, String(power));
            Add(inputs, "+");
        od;
        Add(inputs, "0;\n");
    fi;
    inputs := Concatenation(inputs);
    # Print(inputs);
    # Print(Concatenation(inputs, ReadAll(paricode)));
    # return 0;


    # execute PARI/GP
    Info(InfoAlnuth, 1, "executing PARI/GP with norm_equation.gp");
    output := "";
    Process(DirectoryCurrent(), AL_EXECUTABLE, 
            InputTextString(Concatenation(inputs, ReadAll(paricode))),
            OutputTextString(output,false), 
            Concatenation(["-f", "-q"], [AL_STACKSIZE])
           );

    # close open input stream from file with GP code
    CloseStream(paricode);
    # Placeholder strings:
    x := "E(4)";
    y := ReplacedString("(E(#) + ComplexConjugate(E(#)))", "#", String(cyclotomic_degree));
    # Now parse the output so that GAP can evaluate the square root:
    parsed_output := ReplacedString(output, "x", x);
    parsed_output := ReplacedString(parsed_output, "y", y);
    Print(output, "\n");
    Print(parsed_output, "\n");
    return EvalString(parsed_output);
end);
