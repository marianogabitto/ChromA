from subprocess import Popen, PIPE, call
import multiprocessing
import subprocess
import anndata
import numpy as np
from scipy.sparse import coo_matrix, hstack as sp_hstack, csr_matrix
from scipy.io import mmwrite
import pysam
import copy
import ray
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
            print("Read {} Fragment Files".format(l_))

        assert len(libraries) == len(barcodes)

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
                filename = name + c_.__str__().replace(" ", "").replace("/", "-") + '.tsv'
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


def count_anndata(adata, bed, fragment_slot="fragment_file", barcode_slot="barcodes",
                  unique_id="libraries", name="counts", return_anndata=False):

    print("Parsing Anndata Object for Barcodes, Unique ID & Fragment Files")
    # Validate Presence of Slots in Anndata
    if not ((fragment_slot in adata.obs) and (barcode_slot in adata.obs) and
            (unique_id in adata.obs)):
        print("Fragment or Barcode or Unique ID slot not found in adata.obs")
        return -1
    else:
        # Validate LibraryIDs and
        libraries = adata.obs[fragment_slot].unique()
        library_ids = adata.obs[unique_id].unique()
        if not (len(library_ids) == len(libraries)):
            print("Discordant Libraries Path and IDs")
            return -1

        barcodes = []
        for l_ in libraries:
            working_barcodes = adata.obs[barcode_slot][adata.obs[fragment_slot].isin([l_])].to_numpy()
            barcodes.append(set(working_barcodes))
            print("Read barcodes: {}".format(l_))

        assert len(libraries) == len(barcodes)

    # Parse Bed File/List
    if isinstance(bed, list):
        print("Parsing Bed Object as List of locations")
        if len(bed[0]) < 3:
            print("Bed File contains less than 3 columns")
            return -1
        else:
            interval = bed
            name_bed = "bedlist"
            print("Bed File contains {} regions".format(len(interval)))
    else:
        print("Parsing Bed Object as file")
        if not os.path.exists(bed):
            print("Bed File does not exists")
            return -1
        else:
            interval = list()
            with open(bed, 'r') as f_:
                reader = csv.reader(f_, delimiter='\t')
                for row in reader:
                    if (len(row)) > 0:
                        interval.append([row[0], int(row[1]), int(row[2])])
            name_bed = os.path.split(bed)[1]
            print("Bed File {} contains {} regions".format(name_bed, len(interval)))

    # Setting up containers
    n_cells = np.sum([len(b_) for b_ in barcodes])
    barcode_number = dict()
    count = 0
    for library_barcodes in barcodes:
        for b_ in library_barcodes:
            barcode_number[b_[0]] = count
            count = count + 1
    assert count == n_cells

    # Init Ray
    processors = int(multiprocessing.cpu_count()) - 1
    if not ray.is_initialized():
        ray.init(num_cpus=processors)

    # Split Load by Processor
    splits = int(np.floor(len(interval) / processors))
    split_interval = [interval[i:i + splits] for i in np.arange(0, len(interval), splits, dtype=int)]

    # Split Calculations
    results = []
    for i_ in np.arange(len(split_interval)):
        results.append(counting_many_tsv.remote(libraries, barcode_number, copy.copy(split_interval[i_])))

    # Collect Results
    counts = np.zeros((n_cells, 0), dtype=int)
    # counts = coo_matrix(np.zeros((n_cells, 0), dtype=int))
    for i_ in np.arange(len(split_interval)):
        temp = ray.get(results[i_])
        counts = np.hstack([counts, temp])
        # temp = coo_matrix(ray.get(results[i_]))
        # counts = sp_hstack(counts, temp)

    # #############################################################################
    # Writing Outputs
    print("Formating Outputs")
    counts = csr_matrix(counts)
    print("Outputs")
    #         Writing Count Matrix
    mmwrite(os.path.join(name + "_" + name_bed + "_counts.mtx"), counts)

    #         Writing Barcodes
    f = open(os.path.join(name + "_" + name_bed + "_barcodes.tsv"), "w")
    for i_, library_barcode in enumerate(barcodes):
        for b_ in library_barcode:
            print(b_[0] + "_" + library_ids[i_], end='\n', file=f)
    f.close()

    #         Writing Bed Peaks
    f = open(os.path.join(name + "_" + name_bed + "_peaks.bed"), "w")
    for i, region in enumerate(interval):
        print(region[0], region[1], region[2], end='\n', sep='\t', file=f)
    f.close()

    print("Finished Counting")


@ray.remote(num_cpus=2)
def counting_many_tsv(filenames, barcode_num, inte):
    n_barcodes = np.sum([len(b_) for b_ in barcode_num])
    temp_counts = np.zeros((n_barcodes, len(inte)), dtype=int)

    for i_, f_ in enumerate(filenames):
        if not os.path.exists(f_):
            print("Failed to read bead file")
            return -1
        else:
            sam_file = pysam.TabixFile(f_, threads=2)
        barcodes = barcode_num[i_]

        for j_, b_ in enumerate(inte):
            for read in sam_file.fetch(b_[0], b_[1], b_[2]):
                rows = read.split('\t')
                if barcode_num.__contains__(rows[3]):
                    if (int(rows[1]) > int(b_[1])) and (int(rows[1]) < int(b_[2])):
                        temp_counts[barcodes[rows[3]], j_] += 1
                    if (int(rows[2]) > int(b_[1])) and (int(rows[2]) < int(b_[2])):
                        temp_counts[barcodes[rows[3]], j_] += 1

    return temp_counts
