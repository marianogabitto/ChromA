import multiprocessing
import subprocess
import anndata
import numpy as np
import pysam
import csv
import os


def process_tabix(filename=""):
    if filename == "":
        return

    # Try Create BGZIP
    try:
        out_bgz_fn = filename + ".gz"
        print("Compressing BGZIP File:{}".format(out_bgz_fn))
        with open(out_bgz_fn, 'wb') as out_bgz_fh:
            res = subprocess.call(['bgzip', '-c', filename], stdout=out_bgz_fh)
            if res != 0 or os.stat(out_bgz_fn).st_size == 0:
                print("Could not create BGZIP archive")
                return 0
    except:
        print("Could not create BGZIP archive")
        return 0

    # Try Create TABIX
    try:
        out_index_fn = "{}.tbi".format(out_bgz_fn)
        print("TABIX File:{}".format(out_index_fn))
        if not os.path.exists(out_index_fn):
            subprocess.call(['tabix', out_bgz_fn, '-p', 'bed'])
            if not os.path.exists(out_index_fn) or os.stat(out_index_fn).st_size == 0:
                raise Exception("Error: Could not create index of bgzip archive")
    except:
        print("Could not create TABIX archive")
        return 0


def filter_anndata_barcodes(adata, fragment_slot="fragment_file", barcode_slot="barcodes",
                        write_single_file=False, name="filter", make_tabix=True):
    print("Reading Barcodes & Fragment Files")
    # Validate Presence of Slots in Anndata
    if not ((fragment_slot in adata.obs) and (barcode_slot in adata.obs)):
        print("Fragment or Barcode slot not found in adata.obs")
        return 0
    else:
        libraries = adata.obs[fragment_slot].unique()
        barcodes = []
        for l_ in libraries:
            barcodes.append(adata.obs[barcode_slot][adata.obs[fragment_slot].isin([l_])])

        assert len(libraries) == len(barcodes)
        print("Read {} Fragment Files".format(len(libraries)))

    # FILTER TSV
    filename = ""
    for i_, l_ in enumerate(libraries):
        barc = barcodes[i_]
        if write_single_file:
            filename = name + i_.__str__() + '.tsv'
        else:
            filename = os.path.split(l_)[1] + "_" + name + i_.__str__()  + '.tsv'

        if os.path.splitext(l_)[1] == ".gz":
            print("Filtering File as GZ TABIX: {}".format(l_))
            print("    with {} barcodes".format(len(barc)))

            n_cpu = np.min([16, multiprocessing.cpu_count() - 2])
            tabix_file = pysam.TabixFile(l_, threads=n_cpu)
            with open(filename, 'a') as f_out:
                for j_, row in enumerate(tabix_file.fetch()):
                    print(j_)
                    if row.split('\t')[3] in barc:
                        f_out.write(row + "\n")
            tabix_file.close()

        else:
            print("Filtering File as CSV Equivalent: {}".format(l_))
            print("    with {} barcodes".format(len(barc)))

            with open(l_) as f_in:
                with open(filename, 'a') as f_out:
                    reader = csv.reader(f_in)
                    for row in reader:
                        if row[0].split('\t')[3] in barc:
                            f_out.write(row[0] + "\n")

        if (not write_single_file) and make_tabix:
            process_tabix(filename)

    if write_single_file and make_tabix:
        process_tabix(filename)
