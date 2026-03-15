from SpecialUnitary import *
from sympy.physics.wigner import wigner_6j
from Write_CGCs_Storage import SU2_STORAGE_PATH
import numpy as np
import time

ROOT_PATH = Path(__file__).parent

FLOAT_CUTOFF = 10**(-10)

try:
    import opt_einsum as oe
    disable_opt_einsum = False
except:
    disable_opt_einsum = True

def get_6j_squared_sympy(rep_1 : SU_irrep, rep_2 : SU_irrep, rep_3 : SU_irrep, rep_4 : SU_irrep, rep_5 : SU_irrep, rep_6 : SU_irrep):
    try:
        symbol_6j = wigner_6j(rep_1.get_SU2_j(), rep_2.get_SU2_j(), rep_3.get_SU2_j(), rep_4.get_SU2_j(), rep_5.get_SU2_j(), rep_6.get_SU2_j())
        return float((symbol_6j.evalf())**2)
    except ValueError:
        return 0        # does not fulfill the triangle relation
		

def test_6j_squared_values(young_max_boxes : int):
    """Compare the values obtained for both methods on SU(2): Sympy and CGC method."""

    rep_list = create_representation_list(N=2, horizontal_max=young_max_boxes)

    matches = 0
    errors = 0

    list_size = len(rep_list)
    num_iterations = list_size**6
    progress = 0

    storage = CGC_lists_storage(N=2)
    storage.load_storage(SU2_STORAGE_PATH)

    for rep_1 in rep_list:
        for rep_2 in rep_list:
            for rep_3 in rep_list:
                for rep_4 in rep_list:
                    for rep_5 in rep_list:
                        for rep_6 in rep_list:
                            tensor_valued = get_6j_squared_from_CGC_storage(rep_1, rep_2, rep_3, rep_4, rep_5, rep_6, storage)
                            sympy_valued = get_6j_squared_sympy(rep_1, rep_2, rep_3, rep_4, rep_5, rep_6)
                            
                            #if tensor_valued != 0 and sympy_valued != 0:
                            #    print(f"{tensor_valued} vs {sympy_valued}")#print(abs(tensor_valued - sympy_valued))#print(f"{tensor_valued} vs {sympy_valued}")
                            #    #return

                            if abs(tensor_valued - sympy_valued) < FLOAT_CUTOFF:
                                matches += 1
                            else:
                                errors += 1

                            progress += 1
                        
                        print(f"{(progress / num_iterations):.1%} M:{matches} E:{errors}", end=' \r')
    
    print(f"Matches: {matches}, Errors: {errors}")
	
	

	
def test_6j_squared_performance(young_max_boxes : int):
    """Compare time elapsed for both methods on SU(2): Sympy and CGC method."""

    rep_list = create_representation_list(N=2, horizontal_max=young_max_boxes)
    print(f"Num of reps: {len(rep_list)}")

    list_size = len(rep_list)
    num_iterations = list_size**6
    progress = 0

    # Sympy method

    init_time = time.time()

    for rep_1 in rep_list:
        for rep_2 in rep_list:
            for rep_3 in rep_list:
                for rep_4 in rep_list:
                    for rep_5 in rep_list:
                        for rep_6 in rep_list:
                            #tensor_valued = get_6j_squared(rep_1, rep_2, rep_3, rep_4, rep_5, rep_6)

                            try:
                                symbol_6j = wigner_6j(rep_1.get_SU2_j(), rep_2.get_SU2_j(), rep_3.get_SU2_j(), rep_4.get_SU2_j(), rep_5.get_SU2_j(), rep_6.get_SU2_j())
                                sympy_valued = (symbol_6j.evalf())**2
                            except ValueError:
                                sympy_valued = 0        # does not fulfill the triangle relation

                            progress += 1
                        
                        print(f"{(progress / num_iterations):.1%}", end=' \r')

    final_time = time.time()                        
    
    print(f"Sympy time: {(final_time - init_time):.2f} seconds")

    # Tensor method

    progress = 0

    init_time = time.time()

    storage = CGC_lists_storage(N=2)
    storage.load_storage(SU2_STORAGE_PATH)

    for rep_1 in rep_list:
        for rep_2 in rep_list:
            for rep_3 in rep_list:
                for rep_4 in rep_list:
                    for rep_5 in rep_list:
                        for rep_6 in rep_list:
                            tensor_valued = get_6j_squared_from_CGC_storage(rep_1, rep_2, rep_3, rep_4, rep_5, rep_6, storage)
                            #tensor_valued = get_6j_squared_from_CGCs(rep_1, rep_2, rep_3, rep_4, rep_5, rep_6)
                            progress += 1
                        
                        print(f"{(progress / num_iterations):.1%}", end=' \r')

    final_time = time.time()     

    print(f"Tensor time: {(final_time - init_time):.2f} seconds")



def write_6j_squared(N, young_max_boxes):

    storage = CGC_lists_storage(N)
    storage.load_storage(ROOT_PATH / f"SU{N}_CGC_list.pk")

    symbols_sto = symbols_6j_lists_storage(N, storage)
    symbols_sto.generate_squared_6j_lists(young_max=young_max_boxes, verbose=True)
    symbols_sto.write_storage(ROOT_PATH / f"SU{N}_squared_6j_symbols.pk")


def main():
    write_6j_squared(N=3, young_max_boxes=2)
    #test_6j_squared_values(young_max_boxes=4)
    #test_6j_squared_performance(young_max_boxes=2)




if __name__ == "__main__":
    main()