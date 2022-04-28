from subprocess import Popen, PIPE, call
import multiprocessing
import subprocess
import anndata
import numpy as np
import pysam
import csv
import os


def process_tabix(filename=""):
    """
    Input: tab-separated file, Output: sorted, compressed and tabix file

    """

    # Validate Existance of File and 3-column tab-separated
    if filename == "":
        return
    else:
        with open(filename, 'r') as f_:
            read = csv.reader(f_)
            for row in read:
                if row[0][0] == '#':
                    pass
                else:
                    if len(row[0].split('\t')) < 3:
                        print("Input File cannot be read as a 3-column tab separated file.")
                    else:
                        break

    # Try Create BGZIP and Sorting File
    try:
        out_bgz_fn = filename + ".gz"
        print("Compressing BGZIP File:{}".format(out_bgz_fn))
        with open(out_bgz_fn, 'wb') as out_bgz_fh:
            sort = Popen(["sort", "-k1V", "-k2n", "-k3n", filename], stdout=PIPE)
            bgzip = Popen("bgzip", stdin=sort.stdout, stdout=out_bgz_fh)
            sort.stdout.close()
            _, err = bgzip.communicate()
            if err != None or os.stat(out_bgz_fn).st_size == 0:
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
            call(["tabix", out_bgz_fn, "-s", "1", "-b", "2", "-e", "3"])
            if not os.path.exists(out_index_fn) or os.stat(out_index_fn).st_size == 0:
                raise Exception("Error: Could not create index of bgzip archive")
    except:
        print("Could not create TABIX archive")
        return 0


def filter_anndata_barcodes(adata, fragment_slot="fragment_file", barcode_slot="barcodes",
                            write_single_file=False, name="filter", make_tabix=False, groupby=None):

    print("Reading Barcodes & Fragment Files")

    # Validate Presence of groupby
    if groupby is None:
        print("No Grouping")
        categories = ""
    else:
        if not (groupby in adata.obs):
            print("GroupBy slot not found in adata.obs")
            return 0
        else:
            categories = set(np.unique(adata.obs[groupby]))
            print("Grouping By {} Categories".format(len(categories)))

    # Validate Presence of Slots in Anndata
    if not ((fragment_slot in adata.obs) and (barcode_slot in adata.obs)):
        print("Fragment or Barcode slot not found in adata.obs")
        return 0
    else:
        libraries = adata.obs[fragment_slot].unique()
        barcodes = []
        category_mapping = []
        for l_ in libraries:
            working_barcodes = adata.obs[barcode_slot][adata.obs[fragment_slot].isin([l_])].to_numpy()
            barcodes.append(set(working_barcodes))
            if groupby is None:
                working_categories = np.repeat("", len(working_barcodes))
            else:
                working_categories = adata.obs[groupby][adata.obs[fragment_slot].isin([l_])].to_numpy()
            category_mapping.append(dict(zip(working_barcodes, working_categories)))

        assert len(libraries) == len(barcodes)
        print("Read {} Fragment Files".format(len(libraries)))

    # FILTER TSV
    filename = ""
    for i_, l_ in enumerate(libraries):
        # Access Data from Library
        barc = barcodes[i_]
        mapping = category_mapping[i_]

        # Open Corresponding Output files
        f_out = {}
        for c_ in categories:
            if write_single_file:
                filename = name + i_.__str__() + c_.__str__().replace(" ", "").replace("/", "-") + '.tsv'
            else:
                filename = os.path.split(l_)[1] + "_" + name + \
                           i_.__str__() + c_.__str__().replace(" ", "").replace("/", "-") + '.tsv'
            f_out[c_] = open(filename, 'a+')

        # TODO: TRIVIAL PARALLELIZATION TO SPLIT CALCULATION IN CHROMOSOMES

        if os.path.splitext(l_)[1] == ".gz":
            print("Filtering File as GZ TABIX: {}".format(l_))
            print("    with {} barcodes".format(len(barc)))

            n_cpu = np.min([16, multiprocessing.cpu_count() - 2])
            tabix_file = pysam.TabixFile(l_, threads=n_cpu)

            for j_, row in enumerate(tabix_file.fetch()):
                current_barcode = row.split('\t')[3]
                if current_barcode in barc:
                    category = mapping[current_barcode]
                    f_out[category].write(row + "\n")

            tabix_file.close()

            for k_ in f_out.keys():
                f_out[k_].close()

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
            for c_ in categories:
                filename = os.path.split(l_)[1] + "_" + name + \
                           i_.__str__() + c_.__str__().replace(" ", "").replace("/", "-") + '.tsv'

                print("Making Tabix File of: {}".format(filename))
                process_tabix(filename)

    if write_single_file and make_tabix:
        for c_ in categories:
            filename = name + c_.__str__().replace(" ", "").replace("/", "-") + '.tsv'

            print("Making Single Tabix File of: {}".format(filename))
            process_tabix(filename)
