"""
Computer program based on the paper
arXiv:1009.0437v2 [math-ph] 28 Mar 2011

Author: Valmir Peixôto
"""

import numpy as np
import pickle as pk
from math import sqrt, comb
from pathlib import Path
import copy

flag_warning = False

Q_M_START_INDEX = 0     # it must be 1 as it is defined on the article, but for array purposes it's better to be 0
FLOAT_ZERO_PRECISION = 10**(-10)    # precision to define if a float number is zero

class SU_irrep():
    __slots__ = "i_weight", "N", "dim", "basis", "casimir2", "p_index"

    def __init__(self, irrep_array : list):
        """Initializate a SU(N) irredutible representation."""

        if irrep_array[-1] != 0:
            raise TypeError("Error! Last item of the array must be zero.")

        self.i_weight = irrep_array
        self.N = len(self.i_weight)

        self.dim = None
        self.basis = None
        self.casimir2 = None
        self.p_index = None


    def get_basis(self, crescent_order = True) -> list['SU_state']:
        """Compute the basis for this representation.
        
        It is computed only once. Subsequent calls return the previous result."""
        if self.basis is not None:
            return self.basis
        
        self.basis = create_basis_states_list_for_rep(self, crescent_order=crescent_order)
        return self.basis


    def get_dimension(self) -> int:
        """Calculate the dimension for this SU(N) irreducible representation.
        
        It is computed only once. Subsequent calls return the previous result."""
        if self.dim is not None:
            return self.dim
        
        prod = 1
        for j in range(1, self.N + 1):
            for i in range(1, j):
                frac = ((self.i_weight[i-1] - self.i_weight[j-1])/(j-i))
                prod *= (1 + frac)
            
        self.dim = round(prod)
        return self.dim
    

    def get_quadratic_casimir(self) -> float:
        """Compute the SU(N) quadratic Casimir for this representation.
        
        It is computed only once. Subsequent calls return the previous result."""

        if self.casimir2 is not None:
            return self.casimir2

        num_boxes = sum(self.i_weight)

        sum_value = 0

        # U(N) contribution
        for i in range(0, self.N - 1):
            sum_value += self.i_weight[i] * (self.i_weight[i] - 2*(i+1) + self.N + 1)

        # SU(N) contribution, which is zero for U(N)
        sum_value -= (num_boxes**2) / self.N

        self.casimir2 = sum_value
        return sum_value
    
    def get_conjugate_rep(self) -> 'SU_irrep':
        """Returns the correspondent conjugate representation."""
        column_array = []
        rep_copy = copy.deepcopy(self.i_weight)

        # counting boxes on each column
        while rep_copy[0] > 0:
            boxes = 0
            for i in range(self.N):
                if rep_copy[i] > 0:
                    boxes += 1
                    rep_copy[i] -= 1 
            column_array.append(boxes)

        assert len(column_array) == self.i_weight[0]

        # now apply the conjugate
        for i in range(self.i_weight[0]):
            column_array[i] = self.N - column_array[i]

        # reorder the columns in a way it makes sense for a Young Diagram
        new_column = list(reversed(column_array))

        new_diagram = [ 0 for _ in range(self.N) ]

        for num_boxes in new_column:
            for i in range(num_boxes):
                new_diagram[i] += 1

        return SU_irrep(new_diagram)
    
    def get_p_index(self):
        """Compute P(S) for this representation."""
        if self.p_index is not None:
            return self.p_index
        
        sum_var = 0
        for k in range(1, self.N):
            sum_var += comb(self.N - k + self.i_weight[k-1] - 1, self.N - k)
        self.p_index = sum_var
        return self.p_index

    def generate_highest_state(self, crescent_order = True, compute_zweight = False) -> 'SU_state':
        """Generate the highest state (i.e. it's annihilated by any J+ operator)
        for the representation given."""

        return self.get_basis()[-1] if crescent_order else self.get_basis()[0]
    
    def is_trivial_rep(self):
        return self.i_weight[0] == 0
    
    # SU(2) functions
    
    def get_SU2_j(self) -> float:
        """Get SU(2) corresponding 'j' for this representation."""
        if self.N != 2:
            raise TypeError("Error! This is not a SU(2) representation.")
        return self.i_weight[0] / 2
    
    def get_SU2_j_idx(self) -> int:
        """Get SU(2) corresponding 'j' index for this representation."""
        if self.N != 2:
            raise TypeError("Error! This is not a SU(2) representation.")
        return self.i_weight[0]
    
    # other methods
    
    def __str__(self):
        return str(self.i_weight)
    
    def __eq__(self, other_irrep : 'SU_irrep'):
        return self.i_weight == other_irrep.i_weight
    
    def __lt__(self, other_irrep : 'SU_irrep'):
        if self.N != other_irrep.N:
            raise Exception(f"Trying to compare two representations from different groups:\nSU({self.N}) and SU({other_irrep.N}) respectively.")

        for i in range(self.N):
            if self.i_weight[i] < other_irrep.i_weight[i]:
                return True
            elif self.i_weight[i] > other_irrep.i_weight[i]:
                return False
        return False
    
    def __gt__(self, other_irrep : 'SU_irrep'):
        if self.N != other_irrep.N:
            raise Exception(f"Trying to compare two representations from different groups:\nSU({self.N}) and SU({other_irrep.N}) respectively.")

        for i in range(self.N):
            if self.i_weight[i] > other_irrep.i_weight[i]:
                return True
            elif self.i_weight[i] < other_irrep.i_weight[i]:
                return False
        return False



class SU_state():
    __slots__ = "irrep", "gt_pattern", "N", "qm", "sigma", "z_weight", "p_weight"

    def __init__(self, gt_pattern : list, qm : int = None):
        """Initializate SU(N) state."""
        self.irrep = SU_state.assert_row_order_return_irrep(gt_pattern)
        self.gt_pattern = gt_pattern
        self.N = self.irrep.N
        self.qm = qm
        self.sigma = None
        self.z_weight = None
        self.p_weight = None


    def assert_row_order_return_irrep(gt_pattern : list) -> SU_irrep:
        """Ensure that the rows are GT pattern ordered, and return the irrep the state belongs to."""
        length_list = [ len(row) for row in gt_pattern ]
        if sorted(length_list) != length_list:
            raise TypeError(f"Error! GT pattern is not ordered: {gt_pattern}")
        return SU_irrep(gt_pattern[-1])   # gt_pattern[-1] is the i-weight (or equivalently irrep label)

    
    def get_sigma(self) -> list[int]:
        """Get row sum "sigma" array. 
        
        It is computed only once. Subsequent calls return the previous result."""

        if self.sigma is not None:
            return self.sigma
        
        sigma = []
        num_entries = self.N

        for l in range(num_entries):
            sigma.append(sum(self.gt_pattern[l]))

        self.sigma = sigma
        return sigma
        

    def get_z_weight(self) -> list[float]:
        """Compute the z-weight of the state.
        
        It is computed only once. Subsequent calls return the previous result."""

        if self.z_weight is not None:
            return self.z_weight
        
        z_weight = []

        sigma = self.get_sigma()

        z_weight.append( sigma[0] - 0.5*(sigma[1]) )
        for l in range(1, self.N - 1):
            z_weight.append( sigma[l] - 0.5*(sigma[l+1] + sigma[l-1]) )

        self.z_weight = z_weight
        return z_weight
    

    def get_p_weight(self) -> list[float]:
        """Compute the p-weight of the state.
        
        It is computed only once. Subsequent calls return the previous result."""

        if self.p_weight is not None:
            return self.p_weight
        
        sigma = self.get_sigma()

        self.p_weight = [ sigma[l] - sigma[l-1] if l > 0 else sigma[l] for l in range(self.N) ]
        return self.p_weight


    def get_diagonal_gt_pattern(self):
        """Diagonal type m_{k,l} indices."""
        diagonal_list = []
        pattern_reversed = list(reversed(self.gt_pattern))
        for k in range(self.N):
            diagonal = []    
            for l in range(self.N - k):
                diagonal.append(pattern_reversed[l][k])

            # add Nones for l matches the array index
            # reverse order because highest l is on the top
            diagonal_list.append([None for _ in range(k)] + diagonal[::-1])    

        #print(f"GT-Pattern = {self}, Diagonal pattern = {diagonal_list}")
        return diagonal_list


    def get_qm(self):
        """Get Q(M) for this state."""

        if self.qm is not None:
            return self.qm
        
        basis = self.irrep.get_basis()
        return basis.index(self)


    def compute_j_plus_component(self, arr_k, arr_l) -> float:
        
        k = arr_k + 1       # k = true k
        l = arr_l + 1       # l = true l
        pattern = self.gt_pattern

        prod_1 = 1
        for kp in range(1, l+2):
            prod_1 *= pattern[arr_l+1][kp-1] - pattern[arr_l][arr_k] + k - kp

        prod_2 = 1
        for kp in range(1, l):
            prod_2 *= pattern[arr_l-1][kp-1] - pattern[arr_l][arr_k] + k - kp - 1

        prod_3 = 1
        for kp in range(1, l+1):
            if kp != k:
                num = pattern[arr_l][kp-1] - pattern[arr_l][arr_k] + k - kp
                prod_3 *= num*(num-1)

        assert prod_3 != 0

        return sqrt(-prod_1*prod_2/prod_3)


    def compute_j_minus_component(self, arr_k, arr_l) -> float:
        
        k = arr_k + 1
        l = arr_l + 1
        pattern = self.gt_pattern

        prod_1 = 1
        for kp in range(1, l+2):
            prod_1 *= pattern[arr_l+1][kp-1] - pattern[arr_l][arr_k] + k - kp + 1

        prod_2 = 1
        for kp in range(1, l):
            prod_2 *= pattern[arr_l-1][kp-1] - pattern[arr_l][arr_k] + k - kp

        prod_3 = 1
        for kp in range(1, l+1):
            if kp != k:
                num = pattern[arr_l][kp-1] - pattern[arr_l][arr_k] + k - kp
                prod_3 *= (num+1)*num

        assert prod_3 != 0

        return sqrt(-prod_1*prod_2/prod_3)


    def increased_is_valid(self, k, l) -> bool:
        """Verify if the J+ may create a valid GT-pattern for 'k' and 'l'."""
        
        diagonal_pattern = self.get_diagonal_gt_pattern()

        # verify k == 0 to avoid negative index
        return ((k == 0 or diagonal_pattern[k-1][l-1] >= diagonal_pattern[k][l] + 1 ) 
                    and ( diagonal_pattern[k][l+1] >= diagonal_pattern[k][l] + 1))



    def get_gt_pattern_increment(self, k, l) -> 'SU_state':
        """Get the GT-pattern for J+ acting on this state for 'k' and 'l'."""
        new_state = copy.deepcopy(self.gt_pattern)
        new_state[l][k] += 1
        return SU_state(new_state)
    

    def decreased_is_valid(self, k, l) -> bool:
        """Verify if the J- may create a valid GT-pattern for 'k' and 'l'."""

        diagonal_pattern = self.get_diagonal_gt_pattern()
        
        if diagonal_pattern[k][l] - 1 >= 0:
            # verify k == l for avoiding getting a non-existing entry on GT-pattern
            return ((k == l or diagonal_pattern[k][l] - 1 >= diagonal_pattern[k][l-1])
                    and diagonal_pattern[k][l] - 1 >= diagonal_pattern[k+1][l+1] )
        return False
    

    def get_gt_pattern_decrement(self, k, l) -> 'SU_state':
        """Get the GT-pattern for J- acting on this state for 'k' and 'l'."""
        new_state = copy.deepcopy(self.gt_pattern)
        new_state[l][k] -= 1
        return SU_state(new_state)


    def get_SU2_m(self):
        """Get SU(2) corresponding 'm' for this state."""
        assert self.N == 2
        j = self.irrep.get_SU2_j()
        return self.gt_pattern[0][0] - j
        

    def __str__(self):
        """Convert state rows to string."""
        rows = []
        for np_row in self.gt_pattern[::-1]:
            rows.append(np_row)
        return str(rows)
    

    def __eq__(self, other_state):
        assert self.N == other_state.N
        for i in range(self.N):
            if self.gt_pattern[i] != other_state.gt_pattern[i]:
                return False
        return True

    # TODO: to be tested
    # TODO: is it useful?
    def __lt__(self, other_state):
        """Returns if this state is lying before 'other_state'."""
        N = self.irrep.N

        assert N == other_state.irrep.N

        self_weight = self.get_z_weight()
        othr_weight = other_state.get_z_weight()

        for i in range(N):
            if self_weight[i] - self_weight[N-i] != othr_weight[i] - othr_weight[N-i]:
                return self_weight[i] - self_weight[N-i] < othr_weight[i] - othr_weight[N-i]
        return False
    

class SU_decomposition():
    __slots__ = "decomposition"

    def __init__(self, 
                 rep_1 : SU_irrep, 
                 rep_2 : SU_irrep,
                 decrescent_order = True,
                 rep_aux_list : list[SU_irrep] = None):

        """Decompose two SU(N) representations.

        Arguments:

        'rep_1' : i-weight of irrep 1.
        'rep_2' : i-weight of irrep 2.
        'decrescent_order' : If true, the irreps with highest weights appear first. (default: True)

        The result is stored on the attribute 'decomposition' in a form of a list of 
        tuples (representation, multiplicity).
        """

        if rep_1.N != rep_2.N:
            raise Exception(f"Representations of different groups were given: first one is \na SU({rep_1.N}) representation whereas the second one is a SU({rep_2.N}) representation.")

        initial_array = copy.deepcopy(rep_2.i_weight)
        basis = rep_1.get_basis()

        unique_decomposed_reps_list = []
        multiplicity_list = []

        for state in basis:

            t_array = initial_array.copy()
            decomp_rep = extract_decomposition(state, t_array)     # extract one representation for each state

            # if the state wasn't discarded
            if decomp_rep is not None:
                
                # now verify if it's already on the list
                bFound = False
                for i, previous_rep in enumerate(unique_decomposed_reps_list):
                    
                    # if it's on the list
                    if decomp_rep == previous_rep:

                        multiplicity_list[i] += 1
                        bFound = True
                        break
                
                if not bFound:
                    unique_decomposed_reps_list.append(decomp_rep)
                    multiplicity_list.append(1)

        # the result is on crescent order, so reverse it if needed
        if decrescent_order:
            unique_decomposed_reps_list = unique_decomposed_reps_list[::-1]
            multiplicity_list = multiplicity_list[::-1]

        # init all reps
        unique_decomposed_reps_list = list(map(SU_irrep, unique_decomposed_reps_list))

        if rep_aux_list is not None:
            # considering that "rep_aux_list" have the ordering "top-bottom"
            new_unique_rep_list = []
            max_idx = len(rep_aux_list) - 1
            for unique_rep in unique_decomposed_reps_list:
                p_index = unique_rep.get_p_index()
                if p_index <= max_idx:
                    new_unique_rep_list.append(rep_aux_list[p_index]) # assign a new reference
                else:
                    new_unique_rep_list.append(unique_rep)  # out of list rep

            unique_decomposed_reps_list = new_unique_rep_list

        # now create tuples list
        self.decomposition = list(zip(unique_decomposed_reps_list, multiplicity_list))


    def get_multiplicity(self, rep_final : SU_irrep):
        for decomposed_rep, multiplicity in self.decomposition:
            if rep_final == decomposed_rep:
                return multiplicity
        return 0

    def __str__(self):
        return "\n".join([ f"{rep.i_weight}, multiplicity = {mult}" for rep, mult in self.decomposition ])


def get_fusion_number(rep_1 : SU_irrep,
                          rep_2 : SU_irrep,
                          rep_final : SU_irrep) -> int:
        """Get the multiplicity of 'final_rep' on the decomposition of the tensorial
        product 'rep_1' x 'rep_2'. """
        decomposed_reps_list = SU_decomposition(rep_1, rep_2)
        for decomp_rep, multiplicity in decomposed_reps_list.decomposition:
            if rep_final == decomp_rep:
                return multiplicity
        return 0

class CGC_container(dict):
    """Just a dict that do not create keys for non-existing keys and returns 0.0 instead."""
    def __missing__(self, key):
        return 0.0


class CGC_list():
    __slots__ = "rep_1", "rep_2", "rep_final", "N", "multiplicity", "dim_1", "dim_2", "dim_final", "coefficients"

    def __init__(self, 
                 rep_1 : SU_irrep, 
                 rep_2 : SU_irrep, 
                 rep_final : SU_irrep,
                 multiplicity : int = None,
                 np_dtype = np.double):
        """Create a list of CGCs for the decomposition
        
        rep_1 x rep_2 = multiplicity * rep_final.
        
        If the multiplicity is not furnished, it will be calculated. Otherwise ensure
        giving the correct multiplicity. There are sanity checks for that: errors may be raised
        if the multiplicity is wrong."""

        if not (rep_1.N == rep_2.N and rep_2.N == rep_final.N):        # sanity check
            raise Exception(f"Representations given are from different groups: SU({rep_1.N}), SU({rep_2.N}) and SU({rep_final.N}) respectively.")

        self.coefficients = CGC_container()

        if multiplicity is None:
            multiplicity = get_fusion_number(rep_1, rep_2, rep_final)

        self.rep_1 = rep_1
        self.rep_2 = rep_2
        self.rep_final = rep_final
        self.N = rep_1.N
        self.multiplicity = multiplicity

        self.dim_1 = rep_1.get_dimension()
        self.dim_2 = rep_2.get_dimension()
        self.dim_final = rep_final.get_dimension()

        if multiplicity > 0:

            self.compute_CGC_highest_state()

            #global flag_warning

            #if not flag_warning:
            #    print("Warning: it isnt computing lower CGCs yet!")
            #    flag_warning = True
            #return
        
            # now compute lower states CGCs

            for mult_idx in range(multiplicity):

                done = np.zeros((self.dim_final), dtype=np.int8)
                done[-1] = True     # highest p-weight state is already done

                for qm_final in reversed(range(self.dim_final)):
                    if not done[qm_final]:
                        #print(f"Test : {qm_final}")    # TODO: remove this
                        self.compute_CGC_lower_states(qm_final, mult_idx, done)

    def set_cgc(self,
                qm_1 : int, 
                qm_2 : int, 
                mult_index : int, 
                qm_final : int, 
                value : float):
        """Manually set a CGC coefficient."""
        
        assert 0 <= qm_1 <= self.rep_1.get_dimension()
        assert 0 <= qm_2 <= self.rep_2.get_dimension()
        assert 0 <= mult_index <= self.multiplicity
        assert 0 <= qm_final <= self.rep_final.get_dimension()

        self.coefficients[qm_1, qm_2, mult_index, qm_final] = value


    def compute_CGC_highest_state(self):
        """Computes the list of Clebsch-Gordan Coefficients for the highest weight state on the 
        final representation."""

        def compare_weights(pw_1 : list[float], pw_2 : list[float]):
            """Check if two weights are equal."""
            for i in range(len(pw_1)):
                if abs(pw_1[i] - pw_2[i]) > FLOAT_ZERO_PRECISION:
                    return False
            return True

        if self.multiplicity == 0:
            return 0
        
        dim_1 = self.dim_1
        dim_2 = self.dim_2
        
        # create the tensor basis
        basis_1 = self.rep_1.get_basis()
        basis_2 = self.rep_2.get_basis()

        highest_state = self.rep_final.generate_highest_state()

        dim_final = self.rep_final.get_dimension()

        curr_column = 0     # the column size is the number of CGCs to be computed
        num_states = 0

        # we need to exclude the states which their z-weight sum does not matches the highest state z-weight on final representation
        coeff_mapping = np.full((dim_1, dim_2), -1)
        state_mapping = np.full((dim_1, dim_2), -1)

        HPrimePrime_z_weight = highest_state.get_z_weight()

        #if self.N == 2:
        #    high_m = highest_state.get_SU2_m()
        #    print(f"Highest final rep state: M = {high_m}")    # TODO: SU(2) only
        #else:
        #    high_m = -1
        
        #print(f"Its p-weight: {HPrimePrime_p_weight}")

        #print("####################\n")

        # Setting the number of CGCs to be computed, as well as populate coeff_mapping

        for state_left in basis_1:

            M_z_weight = state_left.get_z_weight()

            for state_right in basis_2:

                MPrime_z_weight = state_right.get_z_weight()

                # verify if there is a possibility of non-vanishing CGCs for these states
                # the sum of z_weight's of the states must matches with the final state z_weight

                sum_z_weight = [ i + j for i, j in zip(M_z_weight, MPrime_z_weight) ]
                


                # TODO: remove this later
                if False:
                    if self.N == 2 and abs(high_m - (state_left.get_SU2_m() + state_right.get_SU2_m())) < FLOAT_ZERO_PRECISION:

                        print(f"Verifying for m1={state_left.get_SU2_m()} and m2={state_right.get_SU2_m()}")  # TODO: SU(2) only
                        print(f"Separate p-weights: {M_z_weight} and {MPrime_z_weight}")
                        print(f"Sum vs Target: {sum_z_weight} vs {HPrimePrime_z_weight}\n")
                    





                if compare_weights(sum_z_weight, HPrimePrime_z_weight):

                    #print(f"Match: {M_z_weight} + {MPrime_z_weight} = {HPrimePrime_z_weight}")

                    coeff_mapping[state_left.qm, state_right.qm] = curr_column
                    curr_column += 1

        #print("####################\n")

        num_cgcs = curr_column
        #print(f"Num of CGCs to be computed: {num_cgcs}")        # TODO: TO BE REMOVED

        assert num_cgcs != 0

        # Ok

        # if it has 1 column (or 1 CGC), then the comparison is 1 to 1
        if num_cgcs == 1:
            for i in range(dim_1):
                for j in range(dim_2):
                    if coeff_mapping[i,j] >= 0:
                        self.set_cgc(i, j, 0, dim_final - 1, 1.0)       # Q(H'') = dim_final - 1
                        return

        # initializate matrix to be solved
        matrix = np.zeros((dim_1 * dim_2, num_cgcs))

        #print(f"Dimension of rep 1: {dim_1}")
        #print(f"Dimension of rep 2: {dim_2}")
        #print(f"Matriz size: {matrix.shape}")
        
        for state_left in basis_1:
            i = state_left.qm

            for state_right in basis_2:
                j = state_right.qm

                if coeff_mapping[i,j] >= 0:
                    
                    # iterate over all J+ operators
                    for l in range(self.N - 1):

                        # iterate over all possible results for J+^l
                        for k in range(l+1):
                            
                            # left J+
                            if state_left.increased_is_valid(k,l):

                                # get corresponding upper state for k and l
                                upper_1 = state_left.get_gt_pattern_increment(k,l)
                                h = basis_1.index(upper_1)  # index on basis 1       # TODO: test performance for .get_qm()

                                if state_mapping[h,j] < 0:
                                    state_mapping[h,j] = num_states
                                    num_states += 1

                                # matrix[row = state][column = cgc]
                                matrix[ state_mapping[h,j] , coeff_mapping[i,j] ] += state_left.compute_j_plus_component(k,l)

                            # right J+
                            if state_right.increased_is_valid(k,l):

                                # get corresponding upper state for k and l
                                upper_2 = state_right.get_gt_pattern_increment(k,l)
                                h = basis_2.index(upper_2)  # index on basis 2       # TODO: test performance for .get_qm()

                                if state_mapping[i,h] < 0:
                                    state_mapping[i,h] = num_states
                                    num_states += 1

                                # matrix[row = state][column = cgc]
                                matrix[ state_mapping[i,h] , coeff_mapping[i,j] ] += state_right.compute_j_plus_component(k,l)

        # matrix is ready
        # solving rectangular matrix by singular value decomposition

        #print(f"matrix to be solved: {matrix}")
        
        u_matrix, singular_values_desc_order, vt_matrix = np.linalg.svd(matrix, compute_uv=True, hermitian=False, full_matrices=True)

        num_zero_singular_values = 0

        # count the zero singular values
        for entry in singular_values_desc_order:
            if abs(entry) < FLOAT_ZERO_PRECISION:
                num_zero_singular_values += 1

        # the number of zero singular values must match the multiplicity, since it has 'multiplicity' linearly independent solutions
        if num_zero_singular_values != self.multiplicity:
            raise Exception(f"The number of zero singular values ({num_zero_singular_values}) does not match the given multiplicity ({self.multiplicity})")

        # now retrieve the CGC's from the singular value decomposition
        for mult_idx in range(self.multiplicity):
            for i in range(dim_1):
                for j in range(dim_2):
                    if coeff_mapping[i,j] >= 0:
                        # get the last rows of V^T (or the last columns of V), which are LI orthonormalized solutions
                        coefficient = vt_matrix[ num_cgcs - mult_idx - 1, coeff_mapping[i,j] ]

                        # verify if it's not zero
                        if abs(coefficient) > FLOAT_ZERO_PRECISION:
                            # the final state is the highest one, so its index is rep_final_dim - 1
                            self.set_cgc(i, j, mult_idx, self.rep_final.get_dimension() - 1, coefficient)


    def compute_CGC_lower_states(self, qm_final, alpha, done):
        """Computes the list of Clebsch-Gordan Coefficients for the lower weight states on the 
        final representation, assuming that the highest weight state CGCs were computed before."""

        state_weight = np.array( self.rep_final.get_basis()[qm_final].get_p_weight() )

        parent_mapping = np.full((self.dim_final), -1)
        multi_mapping = np.full((self.dim_final), -1)
        which_l_mapping = np.full((self.dim_final), -1)

        basis_final = self.rep_final.get_basis()

        num_parents = 0
        num_multi = 0

        for final_rep_state in basis_final:

            p_weight = np.array( final_rep_state.get_p_weight() )

            if np.isclose(p_weight, state_weight).all():    #z_weight == state_weight:
                multi_mapping[final_rep_state.qm] = num_multi
                num_multi += 1
            else:
                #print(p_weight, state_weight)
                for l in range(self.N - 1):  # l = 1, l < N
                    
                    p_weight[l] -= 1
                    p_weight[l+1] += 1

                    if np.isclose(p_weight, state_weight).all(): #z_weight == state_weight:
                        parent_mapping[final_rep_state.qm] = num_parents
                        num_parents += 1
                        which_l_mapping[final_rep_state.qm] = l

                        if not done[final_rep_state.qm]:
                            self.compute_CGC_lower_states(final_rep_state.qm, alpha, done)
                        break

                    p_weight[l] += 1
                    p_weight[l+1] -= 1

        # OBS: (M, N) vs (M, K)

        final_rep_coeffs = np.zeros((num_parents, num_multi))   # coefficients 'b'
        prod_coeffs = np.zeros((num_parents, self.dim_1 * self.dim_2))

        prod_states_mapping = np.full((self.dim_1, self.dim_2), -1)

        num_prod_states = 0

        basis_1 = self.rep_1.get_basis()
        basis_2 = self.rep_2.get_basis()

        for final_rep_state in basis_final:
            if parent_mapping[final_rep_state.qm] >= 0:
                l = which_l_mapping[final_rep_state.qm]

                # left handside of eq.40 of the article
                for k in range(l+1):
                    if final_rep_state.decreased_is_valid(k, l):
                        decreased_state = final_rep_state.get_gt_pattern_decrement(k, l)
                        index = basis_final.index(decreased_state)
                        final_rep_coeffs[ parent_mapping[final_rep_state.qm], multi_mapping[index] ] += final_rep_state.compute_j_minus_component(k, l)

                
                # right handside of eq.40 of the article
                for state_left in basis_1:
                    i = state_left.qm
                    
                    for state_right in basis_2:
                        j = state_right.qm
                        
                        # check for non-vanishing CGC
                        cgc = self.get_CGC(i, j, alpha, final_rep_state.qm)
                        if abs(cgc) > FLOAT_ZERO_PRECISION:
                            for k in range(l+1):
                                
                                # left J+
                                if state_left.decreased_is_valid(k, l):
                                    decreased_state = state_left.get_gt_pattern_decrement(k, l)
                                    index = basis_1.index(decreased_state)

                                    if prod_states_mapping[index,j] < 0:
                                        prod_states_mapping[index,j] = num_prod_states
                                        num_prod_states += 1

                                    prod_coeffs[ parent_mapping[final_rep_state.qm], prod_states_mapping[index,j] ] += cgc * state_left.compute_j_minus_component(k,l)

                                # right J+
                                if state_right.decreased_is_valid(k, l):
                                    decreased_state = state_right.get_gt_pattern_decrement(k, l)
                                    index = basis_2.index(decreased_state)

                                    if prod_states_mapping[i,index] < 0:
                                        prod_states_mapping[i,index] = num_prod_states
                                        num_prod_states += 1

                                    prod_coeffs[ parent_mapping[final_rep_state.qm], prod_states_mapping[i,index] ] += cgc * state_right.compute_j_minus_component(k,l)

        # matrices ready   
        #print(np.zeros((2,2)).shape)
        #print(f"{final_rep_coeffs.shape} vs {prod_coeffs.shape}")
        lstsq_sol, residual, rank, singular_values = np.linalg.lstsq(final_rep_coeffs, prod_coeffs, rcond=None)  # rcond to silence warning
        
        for rep_final_qm in range(self.dim_final):
            if multi_mapping[rep_final_qm] >= 0:
                for i in range(self.dim_1):
                    for j in range(self.dim_2):
                        if prod_states_mapping[i,j] >= 0:
                            cgc = lstsq_sol[ multi_mapping[rep_final_qm], prod_states_mapping[i,j] ]
                            if abs(cgc) > FLOAT_ZERO_PRECISION:
                                self.set_cgc(i,j, alpha, rep_final_qm, cgc)

                done[rep_final_qm] = True

        return


    def get_CGC(self, qm_1, qm_2, mult_index, qm_final):
        if mult_index <= self.multiplicity - 1:
            return self.coefficients[qm_1, qm_2, mult_index, qm_final]
        else:
            raise Exception(f"Requested multiplicity index {mult_index} above maximum: {self.multiplicity - 1}")
    

    def write_CGC(self, filepath : Path | str):
        """Write all CGCs on 'filepath'."""
        if type(filepath) != Path:
            output_path = Path(filepath)
        else:
            output_path = filepath

        with open(output_path, 'wb') as file:
            pk.dump(self.coefficients, file)


    def load_CGC(self, filepath : Path | str):
        """Load CGCs from 'filepath'."""
        if type(filepath) != Path:
            input_path = Path(filepath)
        else:
            input_path = filepath

        with open(input_path, 'rb') as file:
            self.coefficients = pk.load(file)



def sum_decompositions_list(decomp_list : list[tuple[SU_irrep, int]]):

    unique_decomposed_reps_list = []
    multiplicity_list = []

    for decomp in decomp_list:
        for rep, mult in decomp:

            try:
                # suppose it is on the list
                multiplicity_list[ unique_decomposed_reps_list.index(rep) ] += mult
            except ValueError:
                # it is not on the list, so register it
                unique_decomposed_reps_list.append(rep)
                multiplicity_list.append(mult)
                
    return list(zip(unique_decomposed_reps_list, multiplicity_list))


class SU_multiple_decomposition():
    __slots__ = "decomposition"
    def __init__(self, rep_list : list[SU_irrep], decrescent_order = True):

        """Decompose a list of SU(N) representations.

        Arguments:

        'rep_list' : list of the representations to be decomposed.

        The irreps with highest weights always appear first if 'decrescent_order' is true or not given.

        The result is stored on the attribute 'decomposition' in a form of a list of 
        tuples (representation, multiplicity).
        """

        if not len(rep_list) > 1:
            raise Exception("Not enough representations on this list.")

        # old = keep unchanged until sum iteration finishes
        # new = new iteration decomposition list
        current_rep_list_to_decompose_old = [(rep_list[0], 1)]   # first decomposition: first rep, multiplicity 1

        for i in range(len(rep_list) - 1):  # multiplication iteration

            current_rep_list_to_decompose_new = list(zip([], []))   # clean it

            second_operand = rep_list[i+1]    # the next representation at the right on decomposition

            for current_rep, multiplicity in current_rep_list_to_decompose_old:       # sum iteration

                for _ in range(multiplicity):   # account multiplicity: repeat the process
                    curr_decomp = SU_decomposition(current_rep, second_operand)

                    current_rep_list_to_decompose_new = sum_decompositions_list([current_rep_list_to_decompose_new, 
                                                                        curr_decomp.decomposition])
                
            current_rep_list_to_decompose_old = current_rep_list_to_decompose_new   # update

        # finished!
        self.decomposition = sorted(current_rep_list_to_decompose_new, reverse=decrescent_order)

    def __str__(self):
        return "\n".join([ f"{rep.i_weight}, multiplicity = {mult}" for rep, mult in self.decomposition ])



##  Static functions
##

def create_representation_list(N, horizontal_max, ordering = "top-bottom") -> list[SU_irrep]:
    """Creates a list of SU(N) representations using Young Diagrams.

    horizontal_max: The maximum number of boxes to the right on Young diagrams.
    
    ordering: The basis ordering: "left-right" mode first fills the first row with boxes, then 
    the second row and so on. "top-bottom" mode fills Young diagram columns from top to bottom.
    
    Example on SU(3):
    "left-right":
    (0,0,0)
    (1,0,0)
    (2,0,0)
    ...
    (horizontal_max, 0, 0)
    (1, 1, 0)
    ...

    "top-bottom":
    (0,0,0)
    (1,0,0)
    (1,1,0)
    (2,0,0)
    (2,1,0)
    ...
    """

    representation_list = []

    if ordering == "left-right":

        def __reset__(iter: list, j: int):
            """Flat the first 'j' indexes based on 'iter[j]' value."""
            min_value = iter[j]
            while j >= 0:
                iter[j] = min_value
                j -= 1
            return iter

        iter = [ 0 for _ in range(N) ]

        while True:

            representation_list.append(SU_irrep(iter.copy()))

            # if the Young tableau is maximized for the cutoff horizontal_max, stop
            if sum(iter) == horizontal_max*(N-1):
                break

            if iter[0] != horizontal_max:
                iter[0] += 1
            else:

                for j in range(1, N-1):
                    if iter[j] + 1 <= iter[j-1]:
                        iter[j] += 1
                        iter = __reset__(iter, j)
                        break

        return representation_list
    
    elif ordering == "top-bottom":

        iter = [ 0 for _ in range(N) ]

        first_row_boxes = 0

        while True:

            representation_list.append(SU_irrep(iter.copy()))

            if sum(iter) == horizontal_max*(N-1):
                break

            # failed attempt

            curr_idx = N-1 - 1      # 1 for array indexing, 1 for second to last one
            while curr_idx != 0:
                if iter[curr_idx] + 1 <= iter[curr_idx-1]:
                    iter[curr_idx] += 1
                    break
                else:
                    iter[curr_idx] = 0
                curr_idx -= 1
            
            if curr_idx == 0:   # didnt found index, so increment the first row
                iter[0] += 1
                for i in range(1, N):
                    iter[i] = 0


        return representation_list

    else:
        raise Exception(f"Unknown ordering type: {ordering}")

def create_basis_states_list_for_rep(gt_i_weight : SU_irrep | list, crescent_order = True):
    """Given a SU(N) irreducible representation (i-weight) this function computes
    all the basis states which compose the carrier space of the irrep. 
    
    The first state on the list is the lowest-weight state.
    
    'crescent_order' : Returns a basis in a crescent order (default: True)

    'init_zweigths' : Calculate the z-weight for each state and store it. (default: False)
    """

    rep = gt_i_weight if type(gt_i_weight) == SU_irrep else SU_irrep(gt_i_weight)

    basis = [[rep.i_weight]]   # start from the first row 

    max_num_rows = rep.N
    num_rows = 1

    while num_rows != max_num_rows:

        updated_basis = []
        qm = 1

        for state in basis:
        
            last_row = state[-1]

            element_choices_list = []

            # create the list of all possible choices for each element in the next row, obeying the betweeness relation
            # use reversed because the first state of the list will be the highest-weight one
            for i in range(len(last_row)-1):
                element_choices_list.append(list(reversed(range(last_row[i+1], last_row[i] + 1))))

            # now create all possible rows from the list 
            rows_list = []
            length = len(element_choices_list)
            index_list = [0 for _ in range(length)]
            max_index_list = [len(element_choices_list[i]) for i in range(length)]

            finished_rows = False

            while not finished_rows:

                row = [element_choices_list[set_idx][i] for set_idx, i in enumerate(index_list)]
                rows_list.append(row)

                # iterate next index
                for i in reversed(range(len(index_list))):
                    if index_list[i] + 1 < max_index_list[i]:
                        index_list[i] += 1
                        break
                    if i == 0:
                        finished_rows = True
                        break

                    index_list[i] = 0   # reset

            # all possible next rows in 'rows_list'
            for row in rows_list:
                updated_state = copy.deepcopy(state)
                updated_state.append(row)
                updated_basis.append(updated_state)  # add a new state
            
        # end 'for', all old 'basis' elements have been iterated
        basis = updated_basis
        num_rows += 1

    if crescent_order:
        basis = basis[::-1]  # reverse the order

    # now init class for each state

    states_basis = []
    qm = Q_M_START_INDEX
    for state in basis:
        states_basis.append( 
            SU_state(state[::-1], qm)
        )
        qm += 1

    return states_basis

def compute_bkl(state : SU_state):
    """Compute the elements b_{k,l} needed for decomposing representations."""

    gt_label = state.gt_pattern

    size = state.N

    b = np.zeros((size,size), dtype=int)
    
    for l in range(state.N):
        for k in range(l+1):
            if k < l and l >= 1:
                b[k][l] = gt_label[l][k] - gt_label[l-1][k]
            else:
                b[k][l] = gt_label[l][k]
    
    return b


def extract_decomposition(state : SU_state, t_array : list):
    """Compute the representation decomposed from 'gt_pattern' using 't_array'."""

    b = compute_bkl(state)
    num_diagonals = state.N
    assert num_diagonals == len(t_array)

    for k in range(num_diagonals):
        for rev_l in reversed(range(k, num_diagonals)):
            t_array[rev_l] += b[k][rev_l]

            # check inequality
            for i in range(num_diagonals - 1):
                if t_array[i] < t_array[i+1]:
                    return None
    
    # now normalize such that the last element = 0
    last = t_array[-1]

    if last != 0:
        t_array = [ item - last for item in t_array ]

    return t_array


def is_trivial_rep(rep : SU_irrep):
    return rep.i_weight[0] == 0

def generate_SU2_irrep(j : float) -> SU_irrep:
    """Initializate a SU(2) representation for 'j'."""

    state_array = [round(2*j), 0]
    
    return SU_irrep(state_array)

def generate_SU2_state(j, m) -> SU_state:
    """Initializate a SU(2) state for 'j' and 'm'."""
    assert abs(m) <= j      # sanity check

    state_array = [ [round(j+m)],
                   [round(2*j), 0] ]
    
    return SU_state(state_array)

def get_SU2_CGC(j_1, j_2, j_3, m_1, m_2, m_3) -> float:
    """OBS: Very inefficient function, but useful for testing."""
    rep_1 = generate_SU2_irrep(j_1)
    rep_2 = generate_SU2_irrep(j_2)

    rep_final = generate_SU2_irrep(j_3)

    #print(f"Initial reps: {rep_1}, {rep_2}; final rep: {rep_final}")
    #print(f"j1 = {rep_1.get_SU2_j()}, j2 = {rep_2.get_SU2_j()}, J = {rep_final.get_SU2_j()}")

    my_CGC_list = CGC_list(rep_1, rep_2, rep_final)

    state_1 = generate_SU2_state(j_1, m_1)
    state_2 = generate_SU2_state(j_2, m_2)

    final_state = generate_SU2_state(j_3, m_3)

    #print(f"m1 = {state_1.get_SU2_m()}, m2 = {state_2.get_SU2_m()}, M = {final_state.get_SU2_m()}\n")

    return my_CGC_list.get_CGC(state_1.get_qm(), state_2.get_qm(), 0, final_state.get_qm())


def integral_3_matrices(cgc_list : CGC_list, 
                        m_1 : SU_state, 
                        mp_1 : SU_state, 
                        m_2 : SU_state, 
                        mp_2 : SU_state, 
                        m : SU_state, 
                        mp : SU_state):
    """
    Compute the integral of 3 matrices elements from representations in 'cgc_list'.
    """
    sum_var = 0
    for mult_idx in range(cgc_list.multiplicity):
        sum_var += (cgc_list.get_CGC(m_1.get_qm(), m_2.get_qm(), mult_idx, mp.get_qm())
                    * (cgc_list.get_CGC(mp_1.get_qm(), mp_2.get_qm(), mult_idx, m.get_qm())).conjugate() )
    return sum_var / (cgc_list.rep_final.get_dimension())
