from SpecialUnitary import * #SU_irrep, SU_state, SU_decomposition, CGC_list, generate_SU2_irrep
from sympy.physics.wigner import clebsch_gordan
import numpy as np

############### Test functions

def test_old():
    rep = SU_irrep([2,1,0])
    decomp = SU_decomposition(rep, rep)
    print(decomp.get_multiplicity(SU_irrep([4,2,0])))
    print(decomp.get_multiplicity(SU_irrep([2,1,0])))
    print(decomp.get_multiplicity(SU_irrep([1,1,0])))

    return
    #su3_rep_states = create_basis_states_list_for_rep(rep)
    
    #for state in su3_rep_states:
    #    print(state)
    #print(SU_decomposition(rep,rep))

    if False:
        rep_list = create_representation_list(N=4, horizontal_max=5, gt_scheme=True)

        for rep in rep_list:
            hghstate1 = rep.generate_highest_state()
            hghstate2 = rep.generate_highest_state_old()
            if not(hghstate1 == hghstate2):
                print(f"{hghstate1} vs {hghstate2}")
    #else:
        my_rep = SU_irrep([2,1,0])

        basis = my_rep.get_basis() #create_basis_states_list_for_rep()
        for state in basis:
            print(f"{state}, Q(M) = {state.qm}")
        
        #print(my_rep.generate_highest_state())
        #print(basis[-1])

    #rep_list = create_representation_list(N=4, horizontal_max=5)
    

    # SU(2) testing
    j1 = 3/2
    m1 = 3/2

    j2 = 3/2
    m2 = 3/2

    J = 3
    M = 3


    j1_index = int(j1*2)
    j2_index = int(j2*2)
    J_index = int(J*2)

    m1_index = int(j1 + m1)
    m2_index = int(j2 + m2)
    M_index = int(J + M)

    rep = SU_irrep([j1_index,0])
    rep_final = SU_irrep([J_index,0])

    decomp_obj = SU_decomposition(rep, rep)

    for final_rep, mult in decomp_obj.decomposition:
        if final_rep == rep_final:
            rep_final = final_rep
            break
    

    my_list = CGC_list(rep, rep, rep_final, 1)

    print(my_list.get_CGC(m1_index, m2_index, 0, M_index))

    #for state in su3_rep_states:
    #    print(f"State {state}, Wz = {state.get_z_weight()}")
    
    #state = su3_rep_states[-2]
    
    #print(su3_rep_states[-2])
    #print(compute_bkl(su3_rep_states[-2]))


def test_get_particular_CGC():
    print("Testing CGC...")
    rep = SU_irrep([1,0,0])
    decomp = SU_decomposition(rep, rep)

    rep_final = SU_irrep([1,1,0])

    bFound = False
    multiplicity = -1

    #my_state = SU_state([[0], [1,0], [1,0,0]])
    #print(my_state.get_diagonal_gt_pattern())
    #return

    #for decomposed, mult in decomp.decomposition:
    #    print(f"{decomposed}, mult = {mult}")
    #return

    for decomposed, mult in decomp.decomposition:
        if rep_final == decomposed:
            bFound = True
            multiplicity = mult
    
    if not bFound:
        raise Exception("Rep final not found in decomposition")

    my_list = CGC_list(rep, rep, rep_final, multiplicity)

    highest_state_final = rep_final.generate_highest_state()

    state_1 = SU_state([[0], [1,0], [1,0,0]])
    state_2 = SU_state([[1], [1,0], [1,0,0]])

    m1_index = state_1.get_qm()
    m2_index = state_2.get_qm()

    m_index = highest_state_final.get_qm()

    print(my_list.get_CGC(m1_index, m2_index, 0, m_index))


def test_SU2_specific_CGC():
    print("Testing CGC...")

    j1, m1, j2, m2, J, M    =   2, 2    ,    3/2 , 1/2    ,    5/2 , 5/2


    rep_1 = generate_SU2_irrep(j1)
    rep_2 = generate_SU2_irrep(j2)

    rep_final = generate_SU2_irrep(J)

    #print(f"Initial reps: {rep_1}, {rep_2}; final rep: {rep_final}")
    print(f"j1 = {rep_1.get_SU2_j()}, j2 = {rep_2.get_SU2_j()}, J = {rep_final.get_SU2_j()}")

    my_CGC_list = CGC_list(rep_1, rep_2, rep_final)

    state_1 = generate_SU2_state(j1, m1)
    state_2 = generate_SU2_state(j2, m2)

    final_state = generate_SU2_state(J, M)

    print(f"m1 = {state_1.get_SU2_m()}, m2 = {state_2.get_SU2_m()}, M = {final_state.get_SU2_m()}\n")

    my_cgc = my_CGC_list.get_CGC(state_1.get_qm(), state_2.get_qm(), 0, final_state.get_qm())
    print(f"CGC = {my_cgc}, squared = {my_cgc**2}")




def test_SU2_CGCs(j1, j2, J):

    #get_SU2_CGC(2, 3/2, )

    #j1,     j2,     J = (
    #2,      3/2,    5/2)

    rep_1 = generate_SU2_irrep(j1)
    rep_2 = generate_SU2_irrep(j2)

    rep_final = generate_SU2_irrep(J)
    
    print(f"\nCreating CGCs for j1 = {rep_1.get_SU2_j()}, j2 = {rep_2.get_SU2_j()}, J = {rep_final.get_SU2_j()}")

    my_CGC_list = CGC_list(rep_1, rep_2, rep_final)

    basis_1 = rep_1.get_basis()
    basis_2 = rep_2.get_basis()
    basis_final = rep_final.get_basis()

    errors = 0
    matches = 0

    for state_1 in basis_1:
        qm_1 = state_1.get_qm()

        for state_2 in basis_2:
            qm_2 = state_2.get_qm()

            for state_final in basis_final:
                qm_final = state_final.get_qm()

                my_CGC = my_CGC_list.get_CGC(qm_1, qm_2, 0, qm_final)
                if my_CGC != 0:
                    m1 = state_1.get_SU2_m()
                    m2 = state_2.get_SU2_m()
                    M = state_final.get_SU2_m()

                    sympy_cgc = clebsch_gordan(j1, j2, J, m1, m2, M).evalf()

                    if abs(sympy_cgc - my_CGC) > FLOAT_ZERO_PRECISION:   # test includes sign convention
                        print(f"m1 = {m1}, m2 = {m2}, M = {M}, OG: {abs(sympy_cgc):.6f}, My: {abs(my_CGC):.6f}")
                        errors += 1
                    else:
                        matches += 1
    
    print(f"\nMatches: {matches}")
    print(f"Errors: {errors}")






def test_error_assert_exceptions():
    """Catch possible exceptions on CGC calculation..."""
    su2_reps = create_representation_list(N=2, horizontal_max=15)

    errors = 0

    print("SU(2) testing...")
    print(f"Num of reps: {len(su2_reps)}")
    for rep_1 in su2_reps:
        for rep_2 in su2_reps:
            for rep_final in su2_reps:
                try:
                    CGC_list(rep_1, rep_2, rep_final)
                except:
                    print(f"Fail for {rep_1}, {rep_2} and {rep_final}")
                    errors += 1

    reps_list = create_representation_list(N=3, horizontal_max=4)

    print("\nSU(3) testing...")
    print(f"Num of reps: {len(reps_list)}")
    for rep_1 in reps_list:
        for rep_2 in reps_list:
            for rep_final in reps_list:
                try:
                    CGC_list(rep_1, rep_2, rep_final)
                except:
                    print(f"Fail for {rep_1}, {rep_2} and {rep_final}")
                    errors += 1

    if errors == 0:
        print("\nSuccess!\n")


def test_CGCs():
    #test_error_assert_exceptions()
    #test_SU2_specific_CGC()
    #test_SU2_CGCs(2, 3/2, 7/2)
    #test_SU2_CGCs(4/2, 3/2, 5/2)
    #test_SU2_CGCs(6/2, 5/2, 9/2)
    #test_SU2_CGCs(1, 1/2, 1/2)
    #test_SU2_CGCs(5, 3, 8)
    test_SU2_CGCs(10, 4, 10)
    #test_SU2_CGCs(5, 3, 6)
    #test_SU2_CGCs(5, 3, 5)
    #test_SU2_CGCs(5, 3, 4)


def is_CGC_list_orthogonal_type_1(cgc_list : CGC_list) -> bool:

    errors = 0
    
    for alpha in range(cgc_list.multiplicity):  # alpha
        for alpha_t in range(cgc_list.multiplicity):    # alpha tilde
            for Mpp in cgc_list.rep_final.get_basis():      # M''
                for Mpp_t in cgc_list.rep_final.get_basis():    # M'' tilde
                    
                    # expected result from deltas
                    expected_result = 1 if (Mpp == Mpp_t and alpha == alpha_t) else 0
                    
                    # now test combination
                    sum_var = 0
                    for M in cgc_list.rep_1.get_basis():      # M
                        for Mp in cgc_list.rep_2.get_basis():    # M'

                            sum_var += (cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), alpha, Mpp.get_qm()) 
                                        *  (cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), alpha_t, Mpp_t.get_qm())).conjugate()  )

                    if abs(sum_var - expected_result) > FLOAT_ZERO_PRECISION:
                        print(f"{sum_var:.6f} vs {expected_result} for {cgc_list.rep_1} x {cgc_list.rep_2} = {cgc_list.rep_final}; alpha: {alpha}, alpha_t: {alpha_t}, Mpp:{Mpp.get_qm()}, Mpp_t:{Mpp_t.get_qm()}")
                        return False #errors += 1
    
    return (True if not errors else False)



def is_CGC_list_orthogonal_type_2(cgc_list : CGC_list) -> bool:

    errors = 0
    
    for M in cgc_list.rep_1.get_basis():                    # M
        for Mp in cgc_list.rep_2.get_basis():               # M'
            for M_t in cgc_list.rep_1.get_basis():          # M~
                for Mp_t in cgc_list.rep_2.get_basis():     # M'~
                    
                    # expected result from deltas
                    expected_result = 1 if (M == M_t and Mp == Mp_t) else 0
                    
                    # now test combination
                    sum_var = 0
                    for Mpp in cgc_list.rep_final.get_basis():      # M''
                        for alpha in range(cgc_list.multiplicity):    # alpha

                            sum_var += (cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), alpha, Mpp.get_qm()) 
                                        *  (cgc_list.get_CGC(M_t.get_qm(), Mp_t.get_qm(), alpha, Mpp.get_qm())).conjugate()  )

                    if abs(sum_var - expected_result) > FLOAT_ZERO_PRECISION:
                        #print(f"{sum_var:.6f} vs {expected_result}")
                        return False #errors += 1
    
    return (True if not errors else False)


def test_CGCs_orthogonality():
    print("Testing created CGCs orthogonality...")

    N, hor_max = (3, 3) #(2,5) #(3, 2) #(2, 10)

    reps_list = create_representation_list(N=N, horizontal_max=hor_max)

    matches = 0
    errors = 0

    print(f"Testing CGC's orthogonality for SU({N})...")
    print(f"Num of reps: {len(reps_list)}")

    num_iterations = len(reps_list) ** 2
    iterations = 0

    for rep_1 in reps_list:
        for rep_2 in reps_list:
            for rep_final, mult in SU_decomposition(rep_1, rep_2).decomposition:
                my_cgc_list = CGC_list(rep_1, rep_2, rep_final, mult)
                if is_CGC_list_orthogonal_type_1(my_cgc_list):# and is_CGC_list_orthogonal_type_2(my_cgc_list):
                    matches += 1
                else:
                    errors += 1
                
            iterations += 1
            print(f"Matches: {matches}, Errors: {errors}, Progress: {(iterations/num_iterations):.1%}", end=' \r')
    
    print(f"Matches: {matches}, Errors: {errors}")
    print("Finished!")









#############################################
############ test test-functions ############
#############################################


def orthog_type_1(j1, j2, J):
    errors = 0
    for Mpp in np.arange(-J, J+1, 1):          # M''
        for Mpp_t in np.arange(-J, J+1, 1):    # M'' tilde
            
            # expected result from deltas
            expected_result = 1 if (Mpp == Mpp_t) else 0
            
            # now test combination
            sum_var = 0
            for M in np.arange(-j1, j1+1, 1):      # M
                for Mp in np.arange(-j2, j2+1, 1):    # M'

                    sum_var += (clebsch_gordan(j1, j2, J, M, Mp, Mpp) * clebsch_gordan(j1, j2, J, M, Mp, Mpp_t)).evalf()  #(cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), 0, Mpp.get_qm()) 
                               # *  (cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), 0, Mpp_t.get_qm())).conjugate()  )

            if abs(sum_var - expected_result) > FLOAT_ZERO_PRECISION:
                print(f"{sum_var:.6f} vs {expected_result} ")   # for {cgc_list.rep_1} x {cgc_list.rep_2} = {cgc_list.rep_final}; alpha: {alpha}, alpha_t: {alpha_t}, Mpp:{Mpp.get_qm()}, Mpp_t:{Mpp_t.get_qm()}
                return False #errors += 1
    
    return (True if not errors else False)








def orthog_type_2(j1, j2, J):
    errors = 0
    for M in np.arange(-j1, j1+1, 1):                    # M
        for Mp in np.arange(-j2, j2+1, 1):               # M'
            for M_t in np.arange(-j2, j2+1, 1):           # M~
                for Mp_t in np.arange(-j2, j2+1, 1):       # M'~
                    # expected result from deltas
                    expected_result = 1 if (M == M_t and Mp == Mp_t) else 0
                    
                    # now test combination
                    sum_var = 0
                    for Mpp in np.arange(-J, J+1, 1):      # M''
                        sum_var += (clebsch_gordan(j1, j2, J, M, Mp, Mpp) * clebsch_gordan(j1, j2, J, M_t, Mp_t, Mpp)).evalf()  #(cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), 0, Mpp.get_qm()) 
                                    # *  (cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), 0, Mpp_t.get_qm())).conjugate()  )

                    if abs(sum_var - expected_result) > FLOAT_ZERO_PRECISION:
                        #print(f"{sum_var:.6f} vs {expected_result} ")   # for {cgc_list.rep_1} x {cgc_list.rep_2} = {cgc_list.rep_final}; alpha: {alpha}, alpha_t: {alpha_t}, Mpp:{Mpp.get_qm()}, Mpp_t:{Mpp_t.get_qm()}
                        return False #errors += 1
    
    return (True if not errors else False)


    errors = 0
    
    for M in cgc_list.rep_1.get_basis():                    # M
        for Mp in cgc_list.rep_2.get_basis():               # M'
            for M_t in cgc_list.rep_1.get_basis():          # M~
                for Mp_t in cgc_list.rep_2.get_basis():     # M'~
                    
                    # expected result from deltas
                    expected_result = 1 if (M == M_t and Mp == Mp_t) else 0
                    
                    # now test combination
                    sum_var = 0
                    for Mpp in cgc_list.rep_final.get_basis():      # M''
                        for alpha in range(cgc_list.multiplicity):    # alpha

                            sum_var += (cgc_list.get_CGC(M.get_qm(), Mp.get_qm(), alpha, Mpp.get_qm()) 
                                        *  (cgc_list.get_CGC(M_t.get_qm(), Mp_t.get_qm(), alpha, Mpp.get_qm())).conjugate()  )

                    if abs(sum_var - expected_result) > FLOAT_ZERO_PRECISION:
                        print(f"{sum_var:.6f} vs {expected_result}")
                        return False #errors += 1
    
    return (True if not errors else False)





def test_orthogonality_functions():
    print("Testing orthogonality python-coded functions...")

    N, hor_max = (2,5) #(3, 2) #(2, 10)

    reps_list = create_representation_list(N=N, horizontal_max=hor_max)

    matches = 0
    errors = 0

    print(f"Num of reps: {len(reps_list)}")

    num_iterations = len(reps_list) ** 2
    iterations = 0

    for rep_1 in reps_list:
        for rep_2 in reps_list:
            for rep_final, mult in SU_decomposition(rep_1, rep_2).decomposition:
                if orthog_type_1(rep_1.get_SU2_j(), rep_2.get_SU2_j(), rep_final.get_SU2_j()): #is_CGC_list_orthogonal_type_1(my_cgc_list):# and is_CGC_list_orthogonal_type_2(my_cgc_list):
                    matches += 1
                else:
                    errors += 1
                
            iterations += 1
            print(f"Matches: {matches}, Errors: {errors}, Progress: {(iterations/num_iterations):.1%}", end=' \r')
    
    print(f"Matches: {matches}, Errors: {errors}")
    print("Finished!")






if __name__ == "__main__":
    #test_CGCs()
    #test_get_particular_CGC()

    test_CGCs_orthogonality()
    #test_orthogonality_functions()