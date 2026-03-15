from SpecialUnitary import *

ROOT_PATH = Path(__file__).parent

def write_storages():
    su2_storage = CGC_lists_storage(N=2, young_max=5)
    su2_storage.generate_cgc_lists(verbose=True)
    su2_storage.write_storage(ROOT_PATH / "SU2_CGC_list.pk")

def main():
    write_storages()

if __name__ == "__main__":
    main()