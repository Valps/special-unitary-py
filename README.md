# special-unitary-py
Very simple code in Python that creates and handles SU(N) irreducible representations. It computes tensorial product decomposition, Clebsch-Gordan Coefficients, dimensions, quadratic Casimir and basis states of SU(N).

It's based on the article `arXiv:1009.0437 [math-ph]`.

## Main Classes

- `SU_irrep`: Irreducible representation of SU(N). It's labeled by a `N` integer list (irrep weight, or i-weight, in short) in decreasing order, where the last element is always zero. The element `i` represent the number of boxes on row `i` of the representation respective Young diagram.
- `SU_state`: Eigenstate/Eigenvector of `J_z` on SU(N). It's labeled by a Gelfand-Tsetlin pattern.
- `SU_decomposition`: Decomposes the tensorial product of two `SU_irrep` into a sum of representations, each with its own multiplicity.
- `CGC_list`: Generates and stores all Clebsch-Gordan Coefficients (CGCs) arising from the decomposition of the first two `SU_irrep`s into the third `SU_irrep`. It takes into account the multiplicity for SU(N) with N > 2.
- `CGC_lists_storage`: Stores a list of `CGC_list` given `N` and the maximum/cutoff number of boxes on Young diagram to the right. You can write/load its data on disk as a pickle file.
- `SU_multiple_decomposition`: Suitable for decomposing the product of more than two representations. It receives a list of representations whose product is to be decomposed.
- `symbols_6j_lists_storage`: Stores the squared 6j-symbols for SU(N). You can write/load its data on disk as a pickle file

## SU_irrep methods

The methods with a asterisk means that the respective return value is computed only once. Subsequent calls return the previous result, preventing from computing the same value again.

- `get_basis`*: Create the basis for the respective representation. It's a list of `SU_state`s.
- `get_dimension`*: Get the dimension of the carrier space of the representation.
- `get_quadratic_casimir`*: Get the quadratic Casimir eigenvalue of the representation. Normalization follows Tr(T^aT^b) = delta{a,b}, where `T^a` is the `a`-th generator of representation's Lie algebra.
- `get_conjugate_rep`: Get the conjugate representation for the representation. It always returns the same representation on SU(2).
- `get_p_index`*: Get the positive integer index `P` which univocally correspond to a representation of SU(N). It's determined by equation `C6` on the article.
- `generate_highest_state`: Get the highest-weight state of the representation. The highest weight is the one which is always annihiled by a `J^+` operator.
- `is_trivial_rep`: Fast function to check if the representation is the trivial one (i.e. it have no boxes on Young diagram).
- `get_SU2_j`: If the representation's group is SU(2), then you can use this function to get the spin `j` respective to this representation.

## SU_state methods

- `get_z_weight`*: Get the z-weight of the state.
- `get_p_weight`*: Get the p-weight of the state.
- `get_qm`*: Get the `Q(M)` index for the state. It's determined by equation `C7` on the article.
- `get_SU2_m`: If the representation's group is SU(2), then you can use this function to get the quantum number `m` respective to this state.

## SU_decomposition and SU_multiple_decomposition methods

- `get_multiplicity`: Use this function to get the multiplicity of the given final representation on the decomposition.

## CGC_list methods

- `get_CGC`: Get the CGC corresponding to the respective initial and final states `Q(M)`s and for multiplicity index.
- `write_CGC`: Write the CGC list to a file.
- `load_CGC`: Load the CGC list from a file.

## Special Unitary functions

- `create_representation_list`: Given `N` and the maximum number of Young diagram boxes to the right, creates a list with `SU_irrep`s.
- `create_basis_states_list_for_rep`: Given `SU_irrep`, generate a list of `SU_state`s composing a basis for the representation.
- `generate_SU2_irrep`: Create the SU(2) representation respective to given `j`. Note that `j` must be integer or half-integer.
- `generate_SU2_state`: Create the state respective to given `j` and `m`. Note that `j` and `m` must be integer or half-integer.
- `get_SU2_CGC`: It's equivalent to the `sympy.physics.wigner.clebsch_gordan` function of Sympy. It's very inneficient if you're iterating over `j`s or `m`s, in this case use `CGC_list` and `get_CGC` instead.
- `integral_3_matrices`: It computes the integral of 3 matrix elements of SU(N) on Haar measure, where each matrix is associated with its own representation, and the last matrix (the third representation) reads the inverse group element. Contracting all matrix indices gives the multiplicity of the decomposition of the first and second representations into the third.
- `get_6j_squared_from_CGCs`: It computes the 6j symbol squared for given six SU(N) representations.
- `get_6j_squared_from_CGC_storage`: It does the same of the above function but it uses `CGC_lists_storage`. It's recommended to use it if you're iterating over the SU(N) representations, since it's a lot faster than `get_6j_squared_from_CGCs`.
