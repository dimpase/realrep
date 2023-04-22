LoadPackage( "repsn" );
Read("gap_to_pari_norms.g");

DeclareGlobalFunction( "is_real_realisable" );
InstallGlobalFunction( is_real_realisable, function(given_character, the_group)
    local g, sum_fsi, value;
    # Check if the character is real-valued. If it isn't, we don't need to check FSI below, saving time.
    for value in given_character do
        if ImaginaryPart(value) <> 0 then
            return false;
        fi;
    od;
    # If that passes, check its Frobenius-Schur indicator is 1, which is equivalent to this character being real-realisable over this group.
    sum_fsi := 0;
    # -------------- There is an "indicator" function in GAP --------------
    # *** Possible improvement with BSGS? (page 17 K.H. & D.P.) ***
    for g in the_group do
        sum_fsi := sum_fsi + (g*g)^given_character;
    od;
    return sum_fsi = Size(the_group);  # No dividing, just check against the group size.
end);

DeclareGlobalFunction( "nonzero_solution" );
InstallGlobalFunction( nonzero_solution, function(dimension_sq, big_transposed_matrix)
    local chosen_rows, leftover_matrix, maybe_solution, row_index, rhs_vector, vector_to_return;
    # In fact, given that our matrices are guaranteed to have det <> 0, even [1..dim] would suffice.
    for row_index in [1..dimension_sq] do
        # Give minus sign from putting on the RHS:
        rhs_vector := -big_transposed_matrix[row_index];
        # Store remaining rows in a matrix:
        chosen_rows := [1..(row_index - 1)];
        Append(chosen_rows, [(row_index+1)..dimension_sq]);
        leftover_matrix := big_transposed_matrix{chosen_rows};
        maybe_solution := SolutionMat(leftover_matrix, rhs_vector);
        if maybe_solution <> fail then
            # Reconstruct the full solution by putting the 1-entry back.
            vector_to_return := maybe_solution{[1..(row_index - 1)]};
            Add(vector_to_return, 1);
            Append(vector_to_return, maybe_solution{[row_index..(dimension_sq-1)]});
            return vector_to_return;
        fi;
    od;
    # This should never be reached.
    return fail;
end);

DeclareGlobalFunction( "create_matrix_p_sigma_m" );
InstallGlobalFunction( create_matrix_p_sigma_m, function(dimension, the_group, given_representation)
    local collected_tensors, g, gg, mat_gens, matrix_m, matrix_sigma, r, remainder, rep_module, squared_dimension, tensor_minus_identity, vector_form_m, vector_index;
    # Initialize the matrices:
    matrix_m := NullMat(dimension, dimension);
    matrix_sigma := NullMat(dimension, dimension);
    # Then run through the group to compute the sigma.
    # *** Possible improvement with BSGS (page 17 K.H. & D.P.) ***
    for g in the_group do
        r := g^given_representation;
        matrix_sigma := matrix_sigma + TransposedMat(r) * ComplexConjugate(r);  # Mistake before; r*r_con_tr is incorrect.
    od;
    Assert(4, matrix_sigma=ComplexConjugate(TransposedMat(matrix_sigma)));
    # rep_module := GModuleByMats(mat_gens, CyclotomicField(Exponent(the_group)));
    # symm_squares := MTX.InducedActionFactorModule(TensorProductGModule(rep_module, rep_module), WedgeGModule(rep_module));
    # *** In fact, we can just use the Kronecker product, if we are not factoring out the wedge. ***
    mat_gens := GeneratorsOfGroup(Image(given_representation));
    squared_dimension := dimension * dimension;
    collected_tensors := [];
    for g in mat_gens do
        # Need to transpose the matrices, to satisfy g^T * M * g = M, instead of g * M * g^T = M.
        tensor_minus_identity := KroneckerProduct(TransposedMat(g), TransposedMat(g)) - IdentityMat(squared_dimension);
        Append(collected_tensors, tensor_minus_identity);
    od;
    # Now we want to find a nonzero solution in dim^2 variables. So we guess which to fix on 1 and try to solve the rest.
    vector_form_m := nonzero_solution(squared_dimension, TransposedMat(collected_tensors));
    for vector_index in [0..(squared_dimension - 1)] do
        # *** Row-wise or column-wise? Doesn't matter here, as M is symmetric, but it MIGHT matter. ***
        remainder := vector_index mod dimension;
        matrix_m[1 + (vector_index - remainder) / dimension][1 + remainder] := vector_form_m[1 + vector_index];
        # Append(matrix_m, List(vector_form_m{[(1+row_increment)..(dimension + row_increment)]}));
    od;
    Assert(4, matrix_m=TransposedMat(matrix_m));

    for gg in GeneratorsOfGroup(Image(given_representation)) do
        Assert(4, TransposedMat(gg)*matrix_sigma*ComplexConjugate(gg)=matrix_sigma);
        Assert(4, TransposedMat(gg)*matrix_m*gg=matrix_m);
    od;
    return Inverse(matrix_sigma) * matrix_m;
end);

DeclareGlobalFunction( "create_matrix_p_kronecker" );
InstallGlobalFunction( create_matrix_p_kronecker, function(dimension, representation_generators)
    local collected_tensors, g, i, left_end, matrix_p, remainder, right_end, squared_dimension, tensor_minus_identity, vector_form_p;
    squared_dimension := dimension * dimension;
    collected_tensors := [];
    for g in representation_generators do
        # Want to satisfy g_conj * P * g_inv = P. (for every generator of the representation)
        tensor_minus_identity := KroneckerProduct(ComplexConjugate(g), TransposedMat(Inverse(g))) - IdentityMat(squared_dimension);
        Append(collected_tensors, tensor_minus_identity);
    od;
    # Now we want to find a nonzero solution in dim^2 variables. So we guess which to fix on 1 and try to solve the rest.
    vector_form_p := nonzero_solution(squared_dimension, TransposedMat(collected_tensors));
    matrix_p := [];
    left_end := 1;
    right_end := dimension;
    for i in [1..dimension] do
        # Since the Kronecker Product is done row-wise, we have to reconstruct matrix_p row-wise too.
        # For that, we append dimension-sized slices from the vector to the matrix.
        Append(matrix_p, [vector_form_p{[left_end..right_end]}]);
        left_end := left_end + dimension;
        right_end := right_end + dimension;
    od;
    return matrix_p;
end);

DeclareGlobalFunction( "q_conjugate_representation" );
InstallGlobalFunction( q_conjugate_representation, function(rep_dimension, matrix_p, the_group, matrix_generators_of_rep)
    local gen_conjugated, generating_matrix, mat_p_conj, matrix_q, matrix_y, natural_number, new_gen_mats, order_for_Q, order_for_field, primes_in_exp, row, xi_scalar;
    # Select a good order of roots of unity:
    if IsInt(Exponent(the_group) / 4) then
        order_for_Q := 4;  # We have 'i' in our cyclotomic field.
    else
        primes_in_exp := PrimeDivisors(Size(the_group));
        # Here we may want to just take the second element always, but if G === C_p, that would not work.
        if primes_in_exp[1] = 2 then
            order_for_Q := primes_in_exp[2];
        else
            order_for_Q := primes_in_exp[1];
        fi;
    fi;
    # Pick a matrix Q according to the Lemma 3.3.
    # Can be done deterministically.
    mat_p_conj := ComplexConjugate(matrix_p);
    matrix_q := NullMat(rep_dimension, rep_dimension);
    natural_number := 0;
    while Determinant(matrix_q) = 0 do
        xi_scalar := 1 + natural_number * E(order_for_Q);  # The cyclotomic degree could still be reduced, but computing one exponent is often a lot cheaper than a Conductor(flatten(P)).
        matrix_y := matrix_p * xi_scalar;
        matrix_q := ComplexConjugate(matrix_y) + mat_p_conj * matrix_y;
        natural_number := natural_number + 1;  # Derandomising.
    od;
    Assert(4, matrix_p * matrix_q = ComplexConjugate(matrix_q));
    # Now conjugate all generating matrices and
    # Pack it together into another representation format, just like repsn does.
    new_gen_mats := [];
    # Also compute the LCM of conductors, which is the smallest order of cyclotomic field possible for the GModule.
    order_for_field := 1;
    for generating_matrix in matrix_generators_of_rep do
        gen_conjugated := Inverse(matrix_q) * generating_matrix * matrix_q;
        Add(new_gen_mats, gen_conjugated);
        
        for row in gen_conjugated do
            order_for_field := LcmInt( order_for_field, Conductor( row ) );
        od;
        #Display(gen_conjugated);
        Assert(4, gen_conjugated = ComplexConjugate(gen_conjugated));
    od;
    return GroupHomomorphismByImages(the_group, Group(new_gen_mats), GeneratorsOfGroup(the_group), new_gen_mats);
end);

DeclareGlobalFunction( "make_real_representation_NC" );
InstallGlobalFunction( make_real_representation_NC, function(group_G, real_realisable_character)
    local arep, conductor_p, dimension_over_field, gg, gens_group_G, half_dim, i, matrix_p, matrix_rep_generators, norm_mu, number_generators, row_p, norm_root_mu;
    # Print(real_realisable_character, "\n");
    # Compute a real representation of a given group, satisfying a given character.
    # *** ASSUMES that it is real realisable. Perhaps we would want a non-NC version too? ***
    dimension_over_field := real_realisable_character[1];

    # A real-realisable character of dimension 1 is already a real representation; just return.
    if dimension_over_field = 1 then
        number_generators := [];
        gens_group_G := GeneratorsOfGroup(group_G);
        for gg in gens_group_G do
            Add(number_generators, [[gg^real_realisable_character]]);
        od;
        # *** Which is better:  AsGroup([[[1]],[[-1]]])  OR  Group([ [[-1]] ])  ? ***
        return GroupHomomorphismByImages(group_G, AsGroup([[[1]],[[-1]]]), gens_group_G, number_generators);
    fi;
    arep := IrreducibleAffordingRepresentation(real_realisable_character); # We use it later!
    # ^^^^^^^^^^^ Also, need to call it precisely once, because it's non-deterministic.
    
    # Now compute matrix P using the invariant forms, using the Kronecker-product trick, instead of following the Sigma-M procedure of Lemma 3.1:
    matrix_rep_generators := GeneratorsOfGroup(Image(arep));
    #matrix_p := create_matrix_p_sigma_m(dimension_over_field, group_G, arep);
    matrix_p := create_matrix_p_kronecker(dimension_over_field, matrix_rep_generators);
    for gg in matrix_rep_generators do
        Assert(4, matrix_p*gg=ComplexConjugate(gg)*matrix_p);
    od;
    # norm_mu := (matrix_p * ComplexConjugate(matrix_p))[1][1];  # Because it's scaled identity.
    # Instead of doing the entire multiplication, just compute the [1][1] entry:
    norm_mu := 0;
    for i in [1..dimension_over_field] do
        norm_mu := norm_mu + matrix_p[1][i] * ComplexConjugate(matrix_p[i][1]);
    od;
    Assert(4, norm_mu*IdentityMat(dimension_over_field)=matrix_p * ComplexConjugate(matrix_p));
    #Print(norm_mu, "\n");

    #Display(matrix_p);
    #Display(matrix_p * ComplexConjugate(matrix_p));
    half_dim := (dimension_over_field - 1) / 2;
    # Odd dimension case trick:
    if IsInt(half_dim) then
        matrix_p := matrix_p * (norm_mu ^ half_dim / Determinant(matrix_p));
        return q_conjugate_representation(dimension_over_field, matrix_p, group_G, matrix_rep_generators);
    fi;
    # If that isn't the case, PARI/GP helps with the norm-equation.
    # First try to find the norm root in a smaller cyclotomic field:
    conductor_p := 1;
    for row_p in matrix_p do
        conductor_p := LcmInt(conductor_p, Conductor(row_p));
    od;
    #Print("Norms called with Conductor(P):\n");
    norm_root_mu := PariNorms(conductor_p, norm_mu);
    # The norm-root may fail in this case, hence we need the check.
    if norm_root_mu * ComplexConjugate(norm_root_mu) = norm_mu then
        matrix_p := matrix_p / norm_root_mu;
        Assert(4, matrix_p*ComplexConjugate(matrix_p)=IdentityMat(dimension_over_field));
        return q_conjugate_representation(dimension_over_field, matrix_p, group_G, matrix_rep_generators);
    fi;
    # If it fails, we must look into the bigger cyclotomic field (guaranteed success in this case):
    #Print("Norms called with Exponent(G):\n");
    norm_root_mu := PariNorms(Exponent(group_G), norm_mu);
    #Print(norm_root_mu, "\n");
    Assert(1, norm_root_mu * ComplexConjugate(norm_root_mu) = norm_mu, "########################\n\nNorm computation failed!\n\n########################");
    matrix_p := matrix_p / norm_root_mu;
    Assert(4, matrix_p*ComplexConjugate(matrix_p)=IdentityMat(dimension_over_field));
    return q_conjugate_representation(dimension_over_field, matrix_p, group_G, matrix_rep_generators);
end);

# Here are our two entry-points; Either supply a group to compute all real representations in OR
# Provide a representation that you want to get a real version of.

DeclareGlobalFunction( "all_real_representations_of_group" );
InstallGlobalFunction( all_real_representations_of_group, function(group_G)
    local all_real_reps, characters, index, indicators, max_index, repr;
    # Consider all characters of a group, freshly computed.
    # First compute the Frobenius Schur Indicator for all the characters to then compare to:
    indicators := Indicator(CharacterTable(group_G), 2);
    characters := Irr(group_G);
    max_index := Size(indicators);
    all_real_reps := [];
    for index in [1..max_index] do
        if indicators[index] = 1 then
            Add(all_real_reps, make_real_representation_NC(group_G, characters[index]));
        fi;
    od;
    # Not returning the list at the moment, since it's a giant mess to look at in the CMD.
    Print("\n\nHere are all the real representations (up to isomorphism) of this group:\n\n");
    for repr in all_real_reps do
        Print(repr, "\n\n");
    od;
    # return all_real_reps;
end);

DeclareGlobalFunction( "realify_a_representation" );
InstallGlobalFunction( realify_a_representation, function(provided_representation)
    local base_group, conj_representative, the_character;
    # First extract the group from which the representation stems:
    # group_generators = MappingGeneratorsImages(provided_representation)[1];
    base_group := Group(MappingGeneratorsImages(provided_representation)[1]);
    # Then construct the character for this representation:
    the_character := [];
    for conj_representative in ConjugacyClasses(base_group) do
        Add(the_character, TraceMatrix(Representative(conj_representative) ^ provided_representation));
    od;
    # Print(the_character);
    the_character := Character(CharacterTable(base_group), the_character);
    # Print(the_character);
    # *** We may want to do a Frobenius-Schur-Indicator check here... ***
    if is_real_realisable(the_character, base_group) then
        return make_real_representation_NC(base_group, the_character);
    fi;
    Print("Sorry, this representation is not realisable over the real numbers.\n");
end);
