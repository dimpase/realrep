Read("real_representation.g");

DeclareGlobalFunction( "Interesting" );
InstallGlobalFunction( Interesting, function(group_G, real_realisable_character)
    local conductor_p, dimension_over_field, half_dim, matrix_p, norm_mu, row_p, norm_root_mu;
    # Just some text around to make the comparison look more neat.
    Print("Given character:\n");
    Print(real_realisable_character, "\n");
    Print("is realised in RepSn as:\n");
    Print(IrreducibleAffordingRepresentation(real_realisable_character), "\n\n");
    Print("However, it is real realisable, and here is one such representation:\n\n");
    # Compute a real representation of a given group, satisfying a given character.
    # *** ASSUMES that it is real realisable. Perhaps we would want a non-NC version too? ***
    dimension_over_field := real_realisable_character[1];
    # A real-realisable character of dimension 1 is already a real representation; just return.
    if dimension_over_field = 1 then
        return IrreducibleAffordingRepresentation(real_realisable_character);
    fi;
    # Otherwise compute matrix P using the invariant forms, following Lemma 3.1.
    matrix_p := create_matrix_p(dimension_over_field, group_G, IrreducibleAffordingRepresentation(real_realisable_character));
    norm_mu := (matrix_p * ComplexConjugate(matrix_p))[1][1];  # Because it's scaled identity.
    half_dim := (dimension_over_field - 1) / 2;
    # Odd dimension case trick:
    if IsInt(half_dim) then
        matrix_p := matrix_p * (norm_mu ^ half_dim / Determinant(matrix_p));
        return q_conjugate_representation(real_realisable_character, matrix_p, group_G);
    fi;
    # If that isn't the case, the PARI/GP helps with the norm-equation.
    # First try to find the norm root in a smaller cyclotomic field:
    conductor_p := 1;
    for row_p in matrix_p do
        conductor_p := LcmInt(conductor_p, Conductor(row_p));
    od;
    norm_root_mu := PariNorms(conductor_p, norm_mu);
    if norm_root_mu * ComplexConjugate(norm_root_mu) = norm_mu then
        matrix_p := matrix_p / norm_root_mu;
        return q_conjugate_representation(real_realisable_character, matrix_p, group_G);
    fi;
    # If that fails, we must look into the bigger cyclotomic field:
    norm_root_mu := PariNorms(Exponent(group_G), norm_mu);
    Assert(1, norm_root_mu * ComplexConjugate(norm_root_mu) = norm_mu, "########################\n\nNorm computation failed!\n\n########################");
    matrix_p := matrix_p / norm_root_mu;
    return q_conjugate_representation(real_realisable_character, matrix_p, group_G);
end);

# Load this file into GAP, being in the same folder as the "real_representation.g".
# Then just call Interesting(___group___, ___character___);

# Here are some cases for which the character in RepSn induces a complex representation, but is real-realisable:
# SymmetricGroup(4), Character( CharacterTable( SymmetricGroup( [ 1 .. 4 ] ) ), [ 2, 0, 2, -1, 0 ] )
# SymmetricGroup(5), Character( CharacterTable( SymmetricGroup( [ 1 .. 5 ] ) ), [ 4, -2, 0, 1, 1, 0, -1 ] )
# SymmetricGroup(5), Character( CharacterTable( SymmetricGroup( [ 1 .. 5 ] ) ), [ 4, 2, 0, 1, -1, 0, -1 ] )
# DihedralGroup(IsPermGroup, 10), Character( CharacterTable( Group( [ (1,2,3,4,5), (2,5)(3,4) ] ) ), [ 2, 0, E(5)^2+E(5)^3, E(5)+E(5)^4 ] )
# DihedralGroup(IsPermGroup, 36), Character( CharacterTable( Group( [ ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18), ( 2,18)( 3,17)( 4,16)( 5,15)( 6,14)( 7,13)( 8,12)( 9,11) ] ) ), [ 2, 0, 0, -1, -1, 2, -1, -1, 2, -1, -1, 2 ] )
