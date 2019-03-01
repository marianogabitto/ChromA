import seaborn as sns
import numpy as np
import logging
import pickle
import pysam
import copy
import csv
import os

import matplotlib.pyplot as plt
plt.switch_backend("Agg")
from matplotlib.backends.backend_pdf import PdfPages


# ######################################################################
# VALIDATE AND COMPUTE DATA OBJECTS
def mouse_lens():
    lens = {'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116, 'chr5': 151834684,
            'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213, 'chr9': 124595110, 'chr10': 130694993,
            'chr11': 122082543, 'chr12': 120129022, 'chr13': 120421639, 'chr14': 124902244,
            'chr15': 104043685, 'chr16': 98207768, 'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566,
            'chrX': 171031299, 'chrY': 91744698, 'chrM': 16299}

    return lens


def human_lens():
    lens = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
            'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
            'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540,
            'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983,
            'chr20': 63025520, 'chr22': 51304566, 'chr21': 48129895, 'chrX': 155270560, 'chrY': 59373566,
            'chrM': 16571}

    return lens


def build_logger(verbose_mode, filename=None, supress=False):
    # Configure Root Logger
    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    if filename is None:
        handler = logging.StreamHandler()
    else:
        handler = logging.FileHandler(filename)
    formatter = logging.Formatter('%(asctime)s:  %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if verbose_mode == '1':
        logger.setLevel(logging.INFO)
        if not supress:
            logger.info("Running Chrom in Info Mode")
    elif verbose_mode == '2':
        if not supress:
            logger.setLevel(logging.DEBUG)
        logger.debug("Running Chrom in Debug Mode")
    else:
        print("Running Chrom in Silent Mode")
        logger.disabled = True

    # Config Metrics Logger
    logger_metric = logging.getLogger().getChild('metrics')
    logger_metric.propagate = False
    logger_metric.setLevel(logging.INFO)
    if filename is None:
        mhandler = logging.FileHandler('metrics.log')
    else:
        mhandler = logging.FileHandler(filename[:-4] + "_metrics.log")
    mformatter = logging.Formatter('%(asctime)s:  %(message)s')
    mhandler.setFormatter(mformatter)
    logger_metric.addHandler(mhandler)


def validate_inputs(files=None):

    logger = logging.getLogger()
    logger.info("Validating Inputs")
    if files is None:
        logging.error("ERROR: No input files")

    for f_ in files:
        if not os.path.isfile(f_):
            logging.error("ERROR:File does not exists. {}".format(f_))
            raise SystemExit
        if f_[-3:] == 'bam':
            try:
                chr_reads([f_], 'chr1', 3e7, 3e7 + 100)
            except:
                logging.error("ERROR:Could not read file as BAM format. {}".format(f_))
                raise SystemExit
        else:
            try:
                with open(f_, 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    row = next(reader)
                    assert(row[0] in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                      'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                                      'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'])
                    assert(int(row[1]) > 0)
                    assert(int(row[2]) > 0)

                path, name = os.path.split(f_)
                if not os.path.isfile(os.path.join(path, '{}_reads.npy'.format(name))):
                    logging.info("Preprocessing TSV file: {}".format(f_))
                    with open(f_, 'r') as f:
                        max_size = len(f.readlines())

                    with open(f_, 'r') as f:
                        reads_array = np.zeros(max_size, dtype='|S5, int, int, int')
                        reader = csv.reader(f, delimiter='\t')
                        idx = 0
                        for row in reader:
                            reads_array[idx][0] = row[0]
                            reads_array[idx][1] = row[1]
                            reads_array[idx][2] = row[2]
                            reads_array[idx][3] = row[4]
                            idx += 1

                    np.save(os.path.join(path, '{}_reads.npy'.format(name)), reads_array)

            except:
                logging.error("ERROR:Could not read file as TSV format.{}".format(f_))
                raise SystemExit

    logger.info("Inputs Validated")

    return


def regions_th17(filename=None, species='mouse'):
    # ############################################
    # MOUSE REGIONS
    # Tfrc: chr16:32, 580, 000 - 32, 670, 000
    # CTCF: chr8:105, 610, 000 - 105, 705, 000
    # Alas1: chr9:106, 180, 000 - 106, 250, 000
    # Rpl13a: chr7:45, 098, 000 - 45, 160, 000
    # Actb: chr5:142, 880, 000 - 142, 950, 000
    # Stat3: chr11:100, 849, 000 - 100, 945, 000
    # BATF: chr12:85, 666, 000 - 85, 761, 000
    # Fosl2: chr5:32, 095, 000 - 32, 190, 000
    # IRF4: chr13:30, 732, 000 - 30, 825, 000
    # Rorc: chr3:94, 303, 000 - 94, 399, 000
    # #############################################

    # ############################################
    # HUMAN REGIONS
    # Tfrc: chr3:195, 750, 000 - 195, 850, 000
    # Ctcf: chr16:67, 590, 000 - 67, 691, 000
    # Alas1: chr3:52, 215, 600 - 52, 316, 600
    # Rpl13a: chr19:49, 940, 000 - 50, 041, 000
    # Actb:chr7:5, 487, 000 - 5, 588, 300
    # Stat3: chr17:40, 450, 000 - 40, 550, 440
    # Batf: chr14:75, 950, 000 - 76, 050, 000
    # Fosl2: chr2:28, 570, 000 - 28, 670, 000
    # Irf4: chr6:371, 000 - 471, 000
    # Rab7A: chr3:128, 440, 000 - 128, 540, 000
    # #############################################

    logger = logging.getLogger()

    if species == 'mouse':
        regions_list = [['chr16', 32580000, 32670000],
                        ['chr8', 105610000, 105705000],
                        ['chr9', 106180000, 106250000],
                        ['chr7', 45098000, 45160000],
                        ['chr5', 142840000, 142952000],
                        ['chr11', 100849000, 100945000],
                        ['chr12', 85666000, 85761000],
                        ['chr5', 32095000, 32190000],
                        ['chr13', 30732000, 30825000],
                        ['chr3', 94303000, 94399000]]
    elif species == 'human':
        regions_list = [['chr3', 195750000, 195850000],
                        ['chr16', 67590000, 67691000],
                        ['chr3', 52215600, 52316600],
                        ['chr19', 49940000, 50041000],
                        ['chr7', 5487000, 5588300],
                        ['chr17', 40450000, 40550440],
                        ['chr14', 75950000, 76050000],
                        ['chr2', 28570000, 28670000],
                        ['chr6', 371000, 471000],
                        ['chr3', 128440000, 128540000]]
    else:
        print('Wrong species name: mouse/human')
        regions_list = 1
        quit()

    if filename is None:
        quit()

    logger.info("Obtaining Curated Regions")
    obs_vec = []
    length = []
    start_l = []
    chrom_l = []
    for l_ in np.arange(len(regions_list)):

        chrom = regions_list[l_][0]
        start = regions_list[l_][1]
        end = regions_list[l_][2]

        obs_vec.append(chr_reads(filename, chrom, start, end))
        length.append(end - start)
        start_l.append(start)
        chrom_l.append(chrom)

    logger.info("Computing Observations")
    out_data = np.concatenate(obs_vec)
    # out_data = out_data - 5
    # out_data[out_data < 0] = 0
    # idx = out_data[:, 0, 0] > 20
    # out_data[idx, 0, 0] = 20

    return out_data, length, start_l, chrom_l


def regions_chr(filename=None, chromosome=None, species='mouse', blacklisted=True):

    logger = logging.getLogger()
    # Validate Species
    if species == 'mouse':
        chrom_lens = mouse_lens()
    elif species == 'human':
        chrom_lens = human_lens()
    else:
        chrom_lens = []
        logger.error('Wrong species name: mouse/human')

    # Validate Chromosome Name
    if chromosome is None:
        chr_ = []
        logger.error("Forgot chromosome name. Try chrX")
    elif chromosome in chrom_lens.keys():
        chr_ = chromosome
    else:
        chr_ = []
        logger.error("Wrong chromosome name. Try chrX")

    # Compute Coverage of chromosome
    logger.info(chr_ + ": Computing Coverage")
    reads = chr_reads(filename, chr_, 1, chrom_lens[chr_])

    # Compute Chunks
    logger.info(chr_ + ": Parsing")
    chunks = get_chunks(reads, chr_)

    # Obtained reads and populate data structure
    logger.info(chr_ + ": Getting Reads")
    chrom_l, start_l, length, out_data = reads_from_chunks(chunks, reads, chr_)
    out_data = np.concatenate(out_data)

    if blacklisted:
        logger.info(chr_ + ": Removing Blacklisted Regions")
        # Reading List of Blacklisted Regions
        bl_path = os.getcwd()
        bl = read_bed(bl_path + "/ChromA/data/blacklisted/mm10.blacklist.bed")
        """
        bl = []
        for _ in HTSeq.BED_Reader(
                bl_path + "/data/blacklisted/mm10.blacklist.bed"):
            bl.append([_.iv.chrom, _.iv.start, _.iv.end])
        """
        # Removing Blacklisted Regions
        blacklist_reads(out_data, bl, chrom_l, start_l, length)

    logger.info(chr_ + ": Data Collected")

    return out_data, length, start_l, chrom_l


# ######################################################################
# DATA PROCESSING ROUTINES
def chr_reads(files, chrom, start, end, insert_size=False):
    length = int(end - start)
    n_files = len(files)
    out = np.zeros((length, n_files))
    insert_size_calc = []
    number_reads = 0

    for i_, f_ in enumerate(files):
        if f_[-3:] == 'bam':
            sam_file = pysam.AlignmentFile(f_)
            for read in sam_file.fetch(chrom, start, end):
                if not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped \
                        or read.is_duplicate or read.mapping_quality < 30:
                    continue
                else:
                    left_tn5_start = min(read.reference_start, read.next_reference_start) + 4
                    right_tn5_end = left_tn5_start + abs(read.template_length) - 5
                    if insert_size:
                        insert_size_calc.append(np.abs(right_tn5_end - left_tn5_start))
                        number_reads += 1

                    if (left_tn5_start < end - 1) and (left_tn5_start > start + 1):
                        out[left_tn5_start - int(start), i_] += 1

                    if (right_tn5_end < end - 1) and (right_tn5_end > start + 1):
                        out[right_tn5_end - int(start), i_] += 1

        else:
            path, name = os.path.split(f_)

            if os.path.isfile(os.path.join(path, '{}_reads.npy'.format(name))):
                reads_array = np.load(os.path.join(path, '{}_reads.npy'.format(name)))
                chr_idx = reads_array['f0'] == chrom.encode()
                reads_array = reads_array[chr_idx]
                for row in reads_array:
                    left_tn5_start = int(row[1])
                    right_tn5_end = int(row[2])
                    if insert_size:
                        insert_size_calc.append(np.abs(right_tn5_end - left_tn5_start))
                        number_reads += 1

                    if (left_tn5_start < end - 1) and (left_tn5_start > start + 1):
                        out[left_tn5_start - int(start), i_] += 1

                    if (right_tn5_end < end - 1) and (right_tn5_end > start + 1):
                        out[right_tn5_end - int(start), i_] += 1

            else:
                with open(f_, 'r') as tsv_file:
                    reader = csv.reader(tsv_file, delimiter='\t')
                    for row in reader:
                        if row[0] == chrom:
                            left_tn5_start = int(row[1])
                            right_tn5_end = int(row[2])
                            if insert_size:
                                insert_size_calc.append(np.abs(right_tn5_end - left_tn5_start))
                                number_reads += 1

                            if (left_tn5_start < end - 1) and (left_tn5_start > start + 1):
                                out[left_tn5_start - int(start), i_] += 1

                            if (right_tn5_end < end - 1) and (right_tn5_end > start + 1):
                                out[right_tn5_end - int(start), i_] += 1

    if insert_size:
        return out, insert_size_calc, number_reads
    else:
        return out # np.convolve(out[:, 0], np.ones((100, )) / 100, mode='same')[:, None]


def get_chunks(cov, chro, region_size=200000):

    # Collect Stats
    chunk_type = np.zeros(3)

    # Beginning and End of Chromosome
    idxn0 = np.where(cov > 0)[0]
    chr_start = idxn0[0]
    chr_end = idxn0[-1]

    # Bases with 0-Coverage, intervals of 500 and 50 bp.
    # idx0 = np.where(cov[chr_start:chr_end] == 0)[0]
    diff = np.diff(idxn0)
    i_diff500 = np.where(diff > 500)[0]
    i_diff50 = np.where(diff > 50)[0]

    chunks = []
    start = chr_start
    while (chr_end - start) > region_size:
        flag = 0
        # Look for chunks using empty 500 bp
        for i_ in i_diff500[:-1]:
            if ((idxn0[i_] - start) > region_size) and ((idxn0[i_] - start) < 1.5 * region_size):
                chunks.append([start, idxn0[i_]])
                start = idxn0[i_ + 1]
                flag = 1
                chunk_type[0] += 1
                break
        if flag == 1:
            continue
        # Look for chunks using empty 50 bp
        for i_ in i_diff50[:-1]:
            if ((idxn0[i_] - start) > region_size) and ((idxn0[i_] - start) < 1.5 * region_size):
                chunks.append([start, idxn0[i_]])
                start = idxn0[i_ + 1]
                flag = 1
                chunk_type[1] += 1
                break
        # Look for chunks using minimum coverage in 50 bp
        if flag == 0:
            minimal = np.inf
            end = start + region_size
            for r_ in np.arange(start + region_size, start + 1.1 * region_size, 50, dtype=int):
                if np.sum(cov[r_:r_ + 50]) < minimal:
                    end = r_
                    minimal = np.sum(cov[r_:r_ + 50])
            chunks.append([start, end])
            start = end + 50
            chunk_type[2] += 1

    # Last Chunk
    chunks.append([start, chr_end])
    chunks.sort()

    if logging.getLogger().getEffectiveLevel() == 10:
        logger = logging.getLogger()
        avg_length = np.mean(np.array(chunks)[:, 1] - np.array(chunks)[:, 0])
        logger.debug(chro + ": Chunk Avg Length:{0:2.0f}".format(avg_length))
        logger.debug(chro + ": Chunks created using empty 500bp:{0:2.0f}".format(chunk_type[0]))
        logger.debug(chro + ": Chunks created using empty 50bp:{0:2.0f}".format(chunk_type[1]))
        logger.debug(chro + ": Chunks created using minimum 50bp:{0:2.0f}".format(chunk_type[2]))
        data = []
        for i_ in np.arange(len(chunks) - 1):
            data.append([int(chro[3:]), chunks[i_][1], chunks[i_ + 1][0]])
        write_bed("chunks_{}.bed".format(chro), np.array(data))

        temp = np.array(chunks)[:, 1] - np.array(chunks)[:, 0]
        fig1 = plt.figure()
        sns.distplot(temp, hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 0.5},
                     label='Chunk Size ' + chro)
        pp = PdfPages('0_chunk_size' + chro + '.pdf')
        pp.savefig(fig1)
        pp.close()

    return chunks


def reads_from_chunks(chunks, reads, chromo):

    obs_vec = list()
    length = list()
    start_l = list()
    chrom_l = list()

    for c_ in chunks:
        length.append(c_[1] - c_[0])
        start_l.append(c_[0])
        chrom_l.append(chromo)
        obs_vec.append(reads[c_[0]:c_[1]])

    return chrom_l, start_l, length, obs_vec


def count_reads(file, species):

    # Get Logger
    logger = logging.getLogger()

    # Validate Species
    if species == 'mouse':
        chrom_lens = mouse_lens()
    elif species == 'human':
        chrom_lens = human_lens()
    else:
        chrom_lens = []
        logger.error('Wrong species name: mouse/human')
    number_reads = 0

    for c_ in chrom_lens:
        if file[-3:] == 'bam':
            sam_file = pysam.AlignmentFile(file)
            for read in sam_file.fetch(c_, 1, chrom_lens[c_]):
                if not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped \
                        or read.is_duplicate or read.mapping_quality < 30:
                    continue
                else:
                    number_reads += 1

        else:
            logger.error('Accurate Read Count only works for bam files.')

    return number_reads


# ######################################################################
# BEDFILES ROUTINES
def read_bed(bed_file):
    interval = []

    try:
        with open(bed_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if (len(row)) > 0:
                    row[1] = int(row[1])
                    row[2] = int(row[2])
                    interval.append(row)
    except:
        with open(bed_file, 'r') as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                if (len(row)) > 0:
                    row[1] = int(row[1])
                    row[2] = int(row[2])
                    interval.append(row)

    return interval


def write_bed(filename, data, start=None, end=None, ext=100, merge=500):
    # Format Data
    if start is None and end is None:
        chrom = data[:, 0]
        start = np.int32(data[:, 1]) - ext
        end = np.int32(data[:, 2]) + ext
    else:
        chrom = data

    # Merge Close Peaks
    idx = np.where(np.abs(start[1:] - end[:-1]) < merge)[0]
    if idx.shape[0] > 0:
        for i_ in reversed(idx):
            # TODO: Verify both chromosomes are the same.
            # This is highly unlikely not to occur
            end[i_] = end[i_ + 1]
        end = np.delete(end, idx + 1)
        start = np.delete(start, idx + 1)
        chrom = np.delete(chrom, idx + 1)

    with open(filename, 'w') as f:
        for i in np.arange(chrom.shape[0]):
            f.write("chr" + str(int(chrom[i])) + "\t" + str(int(start[i])) + "\t" + str(int(end[i])) + "\n")


def bed_result(filename, data, start, chrom, threshold=0.5):

    # Initial Parameters
    n_datasets = len(data)
    out_regions = []

    # Loop through datasets
    for l_ in np.arange(n_datasets):
        w_data = data[l_] > threshold
        fst = np.where(w_data & ~ np.insert(w_data, 0, 0)[:-1])[0]
        lst = np.where(w_data & ~ np.append(w_data, 0)[1:])[0]
        fst.append(lst[-1], 0) if lst.shape[0] > fst.shape[0] else 0

        # Format Bed
        num_reg = fst.shape[0]
        reg = np.zeros([num_reg, 3])
        for i_ in np.arange(num_reg):
            reg[i_, :] = [int(chrom[l_]), start[l_] + fst[i_], start[l_] + lst[i_] + 1]
        out_regions.append(reg)

    # write bed
    write_bed(filename, np.concatenate(out_regions))

    return np.concatenate(out_regions)


def bed_result_broad_peaks(filename, data, start, chrom, threshold=0.5):

    if data[0].shape[1] < 2:
        print("Not enough states to annotate bed file.")
        return

    def get_overlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))

    # Initial Parameters
    n_datasets = len(data)
    out_regions = []

    # Loop through datasets
    for l_ in np.arange(n_datasets):
        w_data2 = data[l_][:, 2] > threshold
        fst1 = np.where(w_data2 & ~ np.insert(w_data2, 0, 0)[:-1])[0]
        lst1 = np.where(w_data2 & ~ np.append(w_data2, 0)[1:])[0]
        fst1.append(lst1[-1], 0) if lst1.shape[0] > fst1.shape[0] else 0

        w_data12 = data[l_][:, 1:].sum(axis=1) > threshold
        fst12 = np.where(w_data12 & ~ np.insert(w_data12, 0, 0)[:-1])[0]
        lst12 = np.where(w_data12 & ~ np.append(w_data12, 0)[1:])[0]
        fst12.append(lst12[-1], 0) if lst12.shape[0] > fst12.shape[0] else 0

        # Format Bed
        num_reg = fst12.shape[0]
        num_reg1 = fst1.shape[0]
        for i_ in np.arange(num_reg):
            for r_ in np.arange(num_reg1):
                if get_overlap([fst12[i_], lst12[i_]], [fst1[r_], lst1[r_]]) > 0:
                    out_regions.append([int(chrom[l_]), start[l_] + fst12[i_], start[l_] + lst12[i_] + 1])
                    break

    # write bed
    write_bed(filename, np.array(out_regions))

    return np.array(out_regions)


def blacklist_reads(data, bl, chrom, start, length):

    # Find any bl peak within reads and zero-ing
    cum = np.insert(np.cumsum(length)[:-1], 0, 0)
    for i_ in np.arange(len(chrom)):
        for l_ in np.arange(len(bl)):
            if bl[l_][0] == chrom[i_]:
                if (bl[l_][1] > start[i_]) and (bl[l_][2] < start[i_] + length[i_]):
                    st_bl = cum[i_] + bl[l_][1] - start[i_]
                    end_bl = cum[i_] + bl[l_][2] - start[i_]
                    data[st_bl:end_bl, :] = 0


# ######################################################################
# FILE METRICS
def frip_sn(annot, spec='mouse', file=None):

    # Validate Species
    if spec == 'mouse':
        chrom_lens = mouse_lens()
        prom = "ChromA/data/promoters/prom_mm10_genes.bed"
    else:
        chrom_lens = human_lens()
        prom = "ChromA/data/promoters/prom_hg19_genes.bed"
    approx_coef = chrom_lens['chr1']/np.sum(list(chrom_lens.values()))

    # Validate File
    if file is None:
        return 0, 0
    elif not os.path.isfile(file[0]):
        return 0, 0

    # Collect Reads Chromosome 1
    reads, ins, reads_r = chr_reads(file, 'chr1', 1, chrom_lens['chr1'], insert_size=True)

    # Count Fraction of Tn5 Binding Events in Peaks
    frip_count = 0.
    for i_, c_ in enumerate(annot):
        if c_[0] == 1:
            frip_count += reads[int(c_[1]):int(c_[2])].sum()

    frip_count = frip_count / np.sum(reads)

    # Calculate Signal To Noise
    bed_peaks = read_bed(prom)
    prom_length = np.abs(bed_peaks[0][2] - bed_peaks[0][1])
    read_prom = np.zeros(prom_length)
    for b_ in bed_peaks:
        if b_[0] == 'chr1':
            read_prom += reads[int(b_[1]):int(b_[2])][:, 0]
    idx = read_prom.argmax()
    if idx < 200 or idx > 3800:
        idx = 1000
    stn = read_prom[idx - 100:idx + 100].sum() / (read_prom[:100] + read_prom[-100:]).sum()

    # Calculate Insert Size Distribution
    ins = np.array(ins)
    ins_calc = (np.sum(ins[(ins < 210) * (ins > 190)]) / np.sum(ins[(ins < 80) * (ins > 60)]))

    fig1 = plt.figure()
    plt.hist(ins, 300)
    ax = plt.axis()
    plt.axis([0, 1000, ax[2], ax[3]])
    _, name = os.path.split(file[0])
    pp = PdfPages(name[:-4] + '_insert_size.pdf')
    pp.savefig(fig1)
    pp.close()

    return frip_count, stn, ins_calc, reads_r / approx_coef


def metrics(filename, annotations=None, species=None):

    # Create Logger
    logger = logging.getLogger('metrics')
    logger.info("DATASET METRICS.")

    # Compute Metrics
    for f_ in filename:
        logger.info("File:{}".format(f_))

        # FRIP
        if annotations is None:
            logger.info("Run ChromA to compute FRIP.")
        else:
            frip_calc, sn, ins_size, number_r = frip_sn(annotations, species, [f_])
            logger.info("FRIP: {0:3.3f}".format(frip_calc))
            logger.info("Signal To Noise: {0:3.3f}".format(sn))
            logger.info("Insert Size Metric (Ratio MonoNuc to Nuc-Free): {0:3.3f}".format(ins_size))
            logger.info("Number of Reads (Extrapolated from Chromosome 1): {0:5.0e}".format(number_r))

# ######################################################################
# COUNT TSV FRAGMENTS IN PEAKS
"""
def count_in_peaks(barcodes, fragments, peaks, species):
    if species == "mouse":
        limit = 20
    elif species == "human":
        limit = 23
    else:
        limit = 0

    for i in np.arange(1, limit):
        pass

    print("Reading Whitelisted Cells")
    cell_list = deque()
    with open(cells_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i > 1:
                cell = str.split(row[0], ",")
                if not (cell[8] == 'None'):
                    cell_list.append(cell)

    # Reading Bed Regions
    print("Reading Reads in Bed Regions")
    try:
        interval = list()
        with open(bed_file, 'r') as f:
            print("Bed File Tab Separated")
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if (len(row)) > 0:
                    # row[1] = int(row[1])
                    # row[2] = int(row[2])
                    interval.append(row)
    except:
        interval = list()
        with open(bed_file, 'r') as f:
            print("Bed File Space Separated")
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                if (len(row)) > 0:
                    # row[1] = int(row[1])
                    # row[2] = int(row[2])
                    interval.append(row)

    # Setting up containers
    n_cells = len(cell_list)
    n_regions = len(interval)
    barcode_number = dict()
    for i, barcode in enumerate(cell_list):
        barcode_number[barcode[0]] = i

        # Detect Chromosome
    chrom = tsv_file.split('chr')[1].split('.tsv')[0]
    chrom = 'chr' + chrom
    print("Chromosome:" + chrom)

    # Fast Indexing Bed File
    bed = np.array(interval)
    for i, c in enumerate(interval):
        if c[0] == chrom:
            break
    shift = i
    bed = bed[bed[:, 0]==chrom, :]
    bed1 = np.int32(bed[:, 1])
    bed2 = np.int32(bed[:, 2])

    # Parsing TSV
    with open(tsv_file, 'r') as f:
        max_size = len(f.readlines())
    print('File max size:' + max_size.__str__())
    tsv_chr = np.zeros(max_size, dtype='S5').astype(str)
    tsv_st = np.zeros(max_size, dtype=int)
    tsv_end = np.zeros(max_size, dtype=int)
    tsv_barcode = np.zeros(max_size, dtype='S20').astype(str)
    with open(tsv_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        count = 0
        for row in reader:
            [tsv_chr[count], tsv_st[count], tsv_end[count], tsv_barcode[count]] = row[0:4]
            count += 1
            if (count % 10000) == 0:
                print(count/max_size)

    # Collecting Reads in regions
    counts = np.zeros((n_cells, n_regions))
    for i_ in np.arange(len(bed1)):
        if (count / max_size % 1000) == 0:
            print(i_ / len(bed1))
        begin = tsv_st > bed1[i_]
        end = tsv_st < bed2[i_]
        where = np.where((1 * begin + 1 * end) == 2)[0]
        idx = [barcode_number[x] if x in barcode_number else [] for x in tsv_barcode[where]]
        ll = np.array([(item, count) for item, count in Counter(idx).items()])
        if len(ll) > 0:
            counts[ll[:, 0], shift + i_] += ll[:, 1]

        begin = tsv_end > bed1[i_]
        end = tsv_end < bed2[i_]
        where = np.where((1 * begin + 1 * end) == 2)[0]
        ll = np.array([(item, count) for item, count in Counter(idx).items()])
        if len(ll) > 0:
            counts[ll[:, 0], shift + i_] += ll[:, 1]

   # Gathering Outputs"
    python3 gather_counts.py $beds fragments.tsv
"""
