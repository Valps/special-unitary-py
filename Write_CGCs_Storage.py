from SpecialUnitary import *

ROOT_PATH = Path(__file__).parent
SU2_STORAGE_PATH = ROOT_PATH / "SU2_CGC_list.pk"
SU3_STORAGE_PATH = ROOT_PATH / "SU3_CGC_list.pk"

def write_SU2_storages():
    su2_storage = CGC_lists_storage(N=2)
    su2_storage.generate_cgc_lists(young_max=14, verbose=True)
    su2_storage.write_storage(SU2_STORAGE_PATH)

def write_SU3_storages():
    su3_storage = CGC_lists_storage(N=3)
    su3_storage.generate_cgc_lists(young_max=4, verbose=True)
    su3_storage.write_storage(SU3_STORAGE_PATH)


def load_SU2_storage():
    storage = CGC_lists_storage(N=2)
    storage.load_storage(SU2_STORAGE_PATH)
    
    rep_1 = generate_SU2_irrep(1)
    rep_2 = generate_SU2_irrep(1)
    rep_final = generate_SU2_irrep(1/2)

    print(SU_decomposition(rep_1, rep_2))

    cgc_list = storage.find(rep_1, rep_2, rep_final)
    if cgc_list is None:
        print("Fail!")
    else:
        print("Success!")

def main():
    #write_SU2_storages()
    #write_SU3_storages()
    load_SU2_storage()

if __name__ == "__main__":
    main()