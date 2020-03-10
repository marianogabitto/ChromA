from scipy.sparse import csr_matrix
from scipy.io import mmwrite
from collections import deque
import multiprocessing
import seaborn as sns
import numpy as np
import logging
import pysam
import copy
import ray
import csv
import os

import matplotlib
gui_env = ['Agg', 'TKAgg', 'GTKAgg', 'Qt4Agg', 'WXAgg']
for gui in gui_env:
    try:
        matplotlib.use(gui, warn=False, force=True)
        from matplotlib import pyplot as plt
        ff = plt.figure()
        plt.close('all')
        break
    except:
        continue

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


def fly_lens():
    lens = {'chr2L': 23513712, 'chr2R': 25286936, 'chr3L': 28110227, 'chr3R': 32079331,
            'chr4': 1348131, 'chrM': 19524, 'chrX': 23542271, 'chrY': 3667352}

    return lens


def ciona_lens():
    lens = {'Chr1': 15704606, 'Chr10': 8327059, 'Chr11': 8696443, 'Chr12': 7853492, 'Chr13': 4340075, 'Chr14': 6288713,
            'Chr2': 9509017, 'Chr3': 11222275, 'Chr4': 7647452, 'Chr5': 7499939, 'Chr6': 4872821, 'Chr7': 8973615,
            'Chr8': 8275968, 'Chr9': 8318069, 'UAContig48': 29809, 'UAContig24': 49794, 'UAContig20': 53721,
            'UAContig26': 49264, 'UAContig19': 58221, 'UAContig18': 60964, 'UAContig10': 160731, 'UAContig45': 30295,
            'UAContig27': 48695, 'UAContig25': 49624, 'UAContig37': 39605, 'UAContig12': 97350, 'UAContig1': 659517,
            'UAContig2': 455181, 'UAContig17': 61212, 'UAContig40': 34855, 'UAContig3': 422458, 'UAContig6': 302016,
            'UAContig4': 355599, 'UAContig5': 347316, 'UAContig46': 30287, 'UAContig7': 283060, 'UAContig23': 50243,
            'UAContig11': 108407, 'UAContig8': 249787, 'UAContig14': 75836, 'UAContig9': 209129, 'UAContig15': 72792,
            'UAContig31': 43137, 'UAContig33': 41808, 'UAContig13': 87933, 'UAContig51': 26623, 'UAContig42': 32965,
            'UAContig38': 37619, 'UAContig47': 29906, 'UAContig53': 23448, 'UAContig22': 51644, 'UAContig28': 47875,
            'UAContig16': 71705, 'UAContig21': 52875, 'UAContig29': 45916, 'UAContig50': 29685, 'UAContig30': 44216,
            'UAContig41': 33945, 'UAContig44': 30449, 'UAContig39': 34958, 'UAContig43': 31566, 'UAContig34': 41180,
            'UAContig35': 41037, 'UAContig52': 24109, 'UAContig49': 29720, 'UAContig36': 40050, 'UAContig32': 41919}

    return lens


def species_chromosomes(species, file=None):

    if species == 'mouse':
        chrom_length = mouse_lens()
    elif species == 'human':
        chrom_length = human_lens()
    elif species == 'fly':
        chrom_length = fly_lens()
    elif species == 'ciona':
        chrom_length = ciona_lens()
    elif species == 'new':
        if file is None:
            logging.error("File parameter cannot be empty with new species.")
            raise SystemExit
        else:
            chrom_length = dict()
            with open(file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    chrom_length[row[0]] = row[1]
    else:
        logging.error("ERROR: Wrong Species. {}".format(species))
        raise SystemExit

    return chrom_length


def species_regions(species):
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
    elif species == 'fly':
        regions_list = [['chr3L', 7329969, 7379929],
                        ['chr3R', 5600000, 5650000],
                        ['chr4', 1050000, 1100000],
                        ['chr4', 670000, 720000],
                        ['chr2L', 10215000, 10265000]]
    elif species == 'ciona':
        regions_list = [['Chr8', 3360000, 3460000]]
    else:
        logging.info("WARNING: Species {} not found.".format(species))
        regions_list = []

    return regions_list


def species_promoters(species):

    chroma_root = os.path.abspath(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
    chrom_lens = species_chromosomes(species)

    if species == 'mouse':
        prom = chroma_root + "/data/promoters/prom_mm10_genes.bed"
        reference_chromosome = 'chr1'
    elif species == 'fly':
        prom = chroma_root + "/data/promoters/prom_dmel6_genes.bed"
        reference_chromosome = 'chr3R'
    elif species == 'human':
        prom = chroma_root + "/data/promoters/prom_hg19_genes.bed"
        reference_chromosome = 'chr1'
    elif species == 'ciona':
        prom = chroma_root + "/data/promoters/prom_ciona_genes.bed"
        reference_chromosome = 'Chr1'
    else:
        prom, reference_chromosome = [], []

    return chrom_lens, prom, reference_chromosome


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
    """
    The purpose of this routine is to validate the existence of the files
    and in certain cases, try to read the correct file format.
    :param files: List of files to be annotated
    :return:
    """

    logger = logging.getLogger()
    logger.info("Validating Inputs")
    if files is None:
        logging.error("ERROR: No input files")

    for f_ in files:
        if not os.path.isfile(f_):
            logging.error("ERROR:File does not exists. {}".format(f_))
            raise SystemExit

        # TRY TO VALIDATE BAM OR TABIX FILE
        if (f_.split('.')[-1] == 'bam') or (f_.split('.')[-1] == 'gz'):
            try:
                pyfile = pysam.AlignmentFile(f_)
                _ = pyfile.header
            except:
                logging.error("ERROR:Could not read file as BAM format. {}".format(f_))
                raise SystemExit

        # TRY TO VALIDATE AS 3 COLUMNS TSV, BED OR BEDGRAPH FILE
        else:
            try:
                with open(f_, 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    row = next(reader)
                    assert (row[0][:3] == 'chr')
                    assert(int(row[1]) > 0)
                    assert(int(row[2]) > 0)

                preprocess_bedlike(f_)

            except:
                logging.error("ERROR:Could not read file as 3 columns TSV/BED/BEDGRAPH format.{}".format(f_))
                raise SystemExit

    logger.info("Inputs Validated")

    return


def validate_chr(filenames, spec, specfile=None, chr_list=None):

    chrom_length = species_chromosomes(spec, specfile)
    if chr_list is None:
        chrom_list = list(chrom_length.keys())
    else:
        chrom_list = chr_list

    chrom_out = copy.copy(chrom_list)

    for chr_ in chrom_list:
        for f_ in filenames:
            # CHECK FILE EXIST
            if not os.path.isfile(f_):
                logging.error("ERROR:File does not exists. {}".format(f_))
                raise SystemExit

            # CHECK NUMBER OF READS IN CHROMOSOME GREATER THAN 100
            try:
                chromosome = 'chr' + str(chr_)
                reads = chr_reads([f_], chromosome, 1, int(chrom_length[chromosome]))
                if np.sum(reads) < 100:
                    chrom_out.remove(chr_)
                    break
            except:
                chrom_out.remove(chr_)
                break

    logger = logging.getLogger()
    logger.info('Running on mouse genome. Chroms:{}'.format(chr_list))

    return chrom_out


def regions_th17(filename=None, species='mouse', dnase=False):

    logger = logging.getLogger()
    if filename is None:
        logger.error('File Not Found.')
        raise SystemExit

    if species is not 'new':
        regions_list = species_regions(species)
        if len(regions_list) > 0:
            logger.info("Obtaining Curated Regions")
            obs_vec = []
            length = []
            start_l = []
            chrom_l = []
            for l_ in np.arange(len(regions_list)):

                chrom = regions_list[l_][0]
                start = regions_list[l_][1]
                end = regions_list[l_][2]

                obs_vec.append(chr_reads(filename, chrom, start, end, dnase=dnase))
                length.append(end - start)
                start_l.append(start)
                chrom_l.append(chrom)

            logger.info("Computing Observations")
            out_data = np.concatenate(obs_vec)
        else:
            out_data, length, start_l, chrom_l = [], [], [], []
    else:
        out_data, length, start_l, chrom_l = [], [], [], []

    return out_data, length, start_l, chrom_l


def regions_chr(filename=None, chromosome=None, species='mouse', specfile=None, blacklisted=True, dnase=False):

    logger = logging.getLogger()
    chrom_lens = species_chromosomes(species, specfile)

    # Validate Chromosome Name
    if chromosome is None:
        logger.error("Forgot chromosome name. Try chrX")
        raise SystemExit
    elif chromosome in chrom_lens.keys():
        chr_ = chromosome
    else:
        logger.error("Wrong chromosome name. Try chrX")
        raise SystemExit

    # Compute Coverage of chromosome
    logger.info(chr_ + ": Computing Coverage")
    reads = chr_reads(filename, chr_, 1, chrom_lens[chr_], dnase=dnase)

    # Compute Chunks
    logger.info(chr_ + ": Parsing")
    chunks = get_chunks(reads, chr_)

    # Obtained reads and populate data structure
    if len(chunks) > 0:
        logger.info(chr_ + ": Getting Reads")
        chrom_l, start_l, length, out_data = reads_from_chunks(chunks, reads, chr_)
        out_data = np.concatenate(out_data)

        if blacklisted:
            logger.info(chr_ + ": Removing Blacklisted Regions")
            # Reading List of Blacklisted Regions
            chroma_root = os.path.abspath(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
            bl = read_bed(chroma_root + "/data/blacklisted/mm10.blacklist.bed")
            # Removing Blacklisted Regions
            blacklist_reads(out_data, bl, chrom_l, start_l, length)

        logger.info(chr_ + ": Data Collected")
    else:
        chrom_l = []
        start_l = []
        length = []
        out_data = []
        logger.info(chr_ + ": No chunks to parse.")

    return out_data, length, start_l, chrom_l


# ######################################################################
# DATA PROCESSING ROUTINES
def chr_reads(files, chrom, start, end, insert_size=False, dnase=False, map_quality=0):
    # Correct Reads for atac assay or dnase
    if dnase is True:
        corr_right = 0
        corr_left = 0
        insert_size = False
    else:
        corr_right = -5
        corr_left = 4

    length = np.int64(end - start)
    n_files = len(files)
    out = np.zeros((length, n_files))
    insert_size_calc = []
    number_reads = 0

    for i_, f_ in enumerate(files):
        # BAM FILE
        if f_[-3:] == 'bam':
            sam_file = pysam.AlignmentFile(f_)
            for read in sam_file.fetch(chrom, start, end):
                if dnase is False:
                    # Assume Reads are Paired-Ended
                    if not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped \
                            or read.is_duplicate or read.mapping_quality < map_quality:
                        continue
                    else:
                        left_tn5_start = min(read.reference_start, read.next_reference_start) + corr_left
                        right_tn5_end = left_tn5_start + abs(read.template_length) + corr_right
                        if insert_size:
                            insert_size_calc.append(np.abs(right_tn5_end - left_tn5_start))
                            number_reads += 1

                        if (left_tn5_start < end - 1) and (left_tn5_start > start + 1):
                            out[left_tn5_start - int(start), i_] += 1

                        if (right_tn5_end < end - 1) and (right_tn5_end > start + 1):
                            out[right_tn5_end - int(start), i_] += 1
                else:
                    # Assume Reads are Single-Ended
                    if read.mapping_quality > map_quality:
                        # Shift Reads to Interval
                        st = np.min([read.reference_start, read.reference_end]) - int(start)
                        nd = np.max([read.reference_start, read.reference_end]) - int(start)

                        # Clip if one end is outside interval
                        st = np.max([st, 0])
                        nd = np.min([nd, end - int(start)])

                        out[st:(st+4), i_] += 1
                        out[(nd-4):nd, i_] += 1
                        """
                        if read.is_reverse:
                            if (read.reference_end < end - 1) and (read.reference_end > start + 1):
                                out[read.reference_end - int(start), i_] += 1
                        else:
                            if (read.reference_start < end - 1) and (read.reference_start > start + 1):
                                out[read.reference_start - int(start), i_] += 1
                        """

        # TABIX FILE
        elif f_[-3:] == '.gz':
            sam_file = pysam.TabixFile(f_)
            for r_ in sam_file.fetch(chrom, start, end):
                read = r_.split('\t')
                left_tn5_start = int(read[1]) + corr_left
                right_tn5_end = int(read[2]) + corr_right

                if insert_size:
                    insert_size_calc.append(np.abs(right_tn5_end - left_tn5_start))
                    number_reads += 1

                if (left_tn5_start < end - 1) and (left_tn5_start > start + 1):
                    out[left_tn5_start - int(start), i_] += 1

                if (right_tn5_end < end - 1) and (right_tn5_end > start + 1):
                    out[right_tn5_end - int(start), i_] += 1

        # TSV/BED FILE
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
        return out  # np.convolve(out[:, 0], np.ones((100, )) / 100, mode='same')[:, None]


def get_chunks(cov, chro, region_size=200000):

    # Collect Stats
    chunk_type = np.zeros(3)

    # Beginning and End of Chromosome
    idxn0 = np.where(cov > 0)[0]
    if len(idxn0) > 0:
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
    else:
        chunks = []
        if logging.getLogger().getEffectiveLevel() == 10:
            logger = logging.getLogger()
            logger.debug(chro + ": no chunks found.")

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
                row[1] = int(row[1])
                row[2] = int(row[2])
                interval.append(row)
    except:
        with open(bed_file, 'r') as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                row[1] = int(row[1])
                row[2] = int(row[2])
                interval.append(row)

    return interval


def write_bed(filename, data, start=None, end=None, ext=100, merge=500, filterpeaks=25):
    """
    Flexible Routine to Write Bed Files. Can Filter by size, Extends, Merge closeby peaks.

    :param filename:
    :param data:
    :param start:
    :param end:
    :param ext:
    :param merge:
    :param filterpeaks:
    :return:
    """

    # Format Input Data
    if start is None and end is None:
        chrom = data[:, 0]
        start = np.int32(data[:, 1]) - ext
        end = np.int32(data[:, 2]) + ext
    else:
        chrom = data

    # Filter Peaks by size
    idx = np.where(np.abs(end - start) < filterpeaks)[0]
    if idx.shape[0] > 0:
        for i_ in reversed(idx):
            end = np.delete(end, i_)
            start = np.delete(start, i_)
            chrom = np.delete(chrom, i_)

    # Merge Closeby Peaks
    idx = np.where(np.abs(start[1:] - end[:-1]) < merge)[0]
    if idx.shape[0] > 0:
        for i_ in reversed(idx):
            # TODO: Verify both chromosomes are the same.
            # This is highly unlikely to occur
            end[i_] = end[i_ + 1]
        end = np.delete(end, idx + 1)
        start = np.delete(start, idx + 1)
        chrom = np.delete(chrom, idx + 1)

    assert (chrom.shape[0] == start.shape[0])
    assert (chrom.shape[0] == end.shape[0])
    assert (start.shape[0] == end.shape[0])

    with open(filename, 'w') as f:
        for i in np.arange(chrom.shape[0]):
            f.write("chr" + str(chrom[i]) + "\t" + str(int(start[i])) + "\t" + str(int(end[i])) + "\n")


def bed_result(filename, data, start, chrom, threshold=0.5, bedext=100, bedmerge=500, filterpeaks=25):
    """
    Converts posterior into a series of intervals and writes them into bed.
    :param filename:
    :param data:
    :param start:
    :param chrom:
    :param threshold:
    :param bedext:
    :param bedmerge:
    :param filterpeaks:
    :return:
    """

    # Initial Parameters
    n_datasets = len(data)
    out_regions = []

    # Loop through datasets
    chr_l = []
    for l_ in np.arange(n_datasets):
        w_data = data[l_] > threshold
        fst = np.where(w_data & ~ np.insert(w_data, 0, 0)[:-1])[0]
        lst = np.where(w_data & ~ np.append(w_data, 0)[1:])[0]
        fst.append(lst[-1], 0) if lst.shape[0] > fst.shape[0] else 0

        # Format Bed
        num_reg = fst.shape[0]
        reg = np.zeros([num_reg, 2])
        for i_ in np.arange(num_reg):
            chr_l.append(chrom[l_])
            reg[i_, :] = [start[l_] + fst[i_], start[l_] + lst[i_] + 1]
        out_regions.append(reg)

    if len(out_regions) > 0:
        reg_out = np.concatenate(out_regions)
        # write bed
        write_bed(filename, data=np.array(chr_l), start=reg_out[:, 0], end=reg_out[:, 1],
                  ext=bedext, merge=bedmerge, filterpeaks=filterpeaks)
        chr_names = np.core.defchararray.add('chr', np.array(chr_l))
        peaks = np.concatenate([chr_names[:, None], np.int64(reg_out)], axis=1)
    else:
        print("No regions to write to Bed File.")
        peaks = []

    return peaks


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


def preprocess_bedlike(file):
    path, name = os.path.split(file)
    if not os.path.isfile(os.path.join(path, '{}_reads.npy'.format(name))):
        logging.info("Preprocessing TSV file: {}".format(file))
        logging.info("Index file as tabix for accelerated performance")
        with open(file, 'r') as f:
            max_size = len(f.readlines())

        with open(file, 'r') as f:
            reads_array = np.zeros(max_size, dtype='|S5, int, int')
            reader = csv.reader(f, delimiter='\t')
            idx = 0
            for row in reader:
                reads_array[idx][0] = row[0]
                reads_array[idx][1] = row[1]
                reads_array[idx][2] = row[2]
                idx += 1

        np.save(os.path.join(path, '{}_reads.npy'.format(name)), reads_array)


# ######################################################################
# FILE METRICS
def frip_sn(annot, spec='mouse', file=None, dnase=False):

    chrom_lens, prom, reference_chromosome = species_promoters(spec)
    approx_coef = chrom_lens[reference_chromosome] / np.sum(list(chrom_lens.values()))

    # Validate File
    if file is None:
        return 0, 0
    elif not os.path.isfile(file[0]):
        return 0, 0

    # Collect Reads Reference Chromosome
    reads, ins, reads_r = chr_reads(file, reference_chromosome, 1,
                                    chrom_lens[reference_chromosome], insert_size=True, dnase=dnase)

    # Count Fraction of Tn5 Binding Events in Peaks
    frip_count = 0.
    for i_, c_ in enumerate(annot):
        if c_[0] == reference_chromosome:
            frip_count += reads[np.int64(c_[1]):np.int64(c_[2])].sum()

    frip_count = frip_count / np.sum(reads)

    # Calculate Signal To Noise
    bed_peaks = read_bed(prom)
    prom_length = np.abs(bed_peaks[0][2] - bed_peaks[0][1])
    read_prom = np.zeros(prom_length)
    for b_ in bed_peaks:
        if b_[0] == reference_chromosome:
            read_prom += reads[int(b_[1]):int(b_[2])][:, 0]
    idx = read_prom.argmax()
    if idx < 200 or idx > 3800:
        idx = 1000
    stn = read_prom[idx - 100:idx + 100].sum() / (read_prom[:100] + read_prom[-100:]).sum()

    # Calculate Insert Size Distribution
    ins = np.array(ins)
    ins_calc = (np.sum(ins[(ins < 210) * (ins > 190)]) / (1 + np.sum(ins[(ins < 80) * (ins > 60)])))

    try:
        fig1 = plt.figure()
        plt.hist(ins, 300)
        ax = plt.axis()
        plt.axis([0, 1000, ax[2], ax[3]])
        name, _ = os.path.splitext(file[0])
        pp = PdfPages(name + '_insert_size.pdf')
        pp.savefig(fig1)
        pp.close()
    except:
        print("Matplotlib ERROR Generating Insert Size Distribution Plot")
        pass

    return frip_count, stn, ins_calc, reads_r / approx_coef


def metrics(filename, annotations=None, species=None):

    # Create Logger
    logger = logging.getLogger('metrics')
    logger.info("DATASET METRICS.")

    # Compute Metrics
    logger.info("DATASET METRICS.")
    if (species is None) or (species is "new"):
        logger.info("Cannot Compute metrics with {} species".format(species))
    else:
        for f_ in filename:
            logger.info("File:{}".format(f_))
            if annotations is None:
                logger.info("Run ChromA to compute metrics.")
            else:
                frip_calc, sn, ins_size, number_r = frip_sn(annotations, species, [f_])
                logger.info("FRIP: {0:3.3f}".format(frip_calc))
                logger.info("Signal To Noise: {0:3.3f}".format(sn))
                logger.info("Insert Size Metric (Ratio MonoNuc to Nuc-Free): {0:3.3f}".format(ins_size))
                logger.info("Number of Reads (Extrapolated from Chromosome 1): {0:5.0e}".format(number_r))


# ######################################################################
# COUNT FRAGMENTS IN PEAKS
def count_fragments_bed(tsv_file, bed_file, cells_file):

    # Reading Correct Cell Barcodes
    print("Reading Whitelisted Cells")
    cell_list = deque()
    filetype = 'cellranger'
    with open(cells_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i == 0:
                header = str.split(row[0], ",")
                # Check if CellRanger Format by length of fields
                if len(header) < 8:
                    filetype = ''
                # If neither 'x' nor 'barcode' is in the first line, it must be a barcode. Save it
                if not(('barcode' in header[0]) or ('x' in header[0])):
                    if filetype == 'cellranger':
                        cell = str.split(row[0], ",")
                        if not (cell[8] == 'None'):
                            cell_list.append(cell)
                    else:
                        cell_list.append(row)

            if i > 0:
                if filetype == 'cellranger':
                    cell = str.split(row[0], ",")
                    if not (cell[8] == 'None'):
                        cell_list.append(cell)
                else:
                    cell_list.append(row)

    # Reading Bed Regions
    print("Reading Reads in Bed Regions")
    try:
        interval = list()
        with open(bed_file, 'r') as f:
            print("tab")
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if (len(row)) > 0:
                    row[1] = int(row[1])
                    row[2] = int(row[2])
                    interval.append(row)
    except:
        interval = list()
        with open(bed_file, 'r') as f:
            print("space ")
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                if (len(row)) > 0:
                    row[1] = int(row[1])
                    row[2] = int(row[2])
                    interval.append(row)

    # Setting up containers
    n_cells = len(cell_list)
    barcode_number = dict()
    for i, barcode in enumerate(cell_list):
        barcode_number[barcode[0]] = i

    # Setting up paths
    path, file = os.path.split(tsv_file)
    path_bed, bfile = os.path.split(bed_file)

    # Init Ray
    processors = int(multiprocessing.cpu_count()) - 1
    if not ray.is_initialized():
        ray.init(num_cpus=processors, object_store_memory=int(40e9), include_webui=False)

    # Split Load by Processor
    splits = int(np.floor(len(interval) / processors))
    split_interval = [interval[i:i + splits] for i in np.arange(0, len(interval), splits, dtype=int)]

    # Split Calculations
    results = []
    for i_ in np.arange(len(split_interval)):
        results.append(filtering_tsv.remote(tsv_file, barcode_number, copy.copy(split_interval[i_])))
        # filtering_tsv(tsv_file, barcode_number, copy.copy(split_interval[i_]))

    # Collect Results
    counts = np.zeros((n_cells, 0), dtype=int)
    for i_ in np.arange(len(split_interval)):
        temp = ray.get(results[i_])
        counts = np.hstack([counts, temp])

    # #############################################################################
    # Writing Outputs
    print("Formating Outputs")
    counts = csr_matrix(counts)
    print("Outputs")
    mmwrite(os.path.join(path, file[:-3] + "_" + bfile[:-4] + "_counts.mtx"), counts)

    #         Writing Barcodes
    f = open(os.path.join(path, file[:-3] + "_" + bfile[:-4] + "_barcodes.tsv"), "w")
    for i, barcode in enumerate(cell_list):
        print(barcode[0], end='\n', file=f)
    f.close()

    #         Writing Bed Peaks
    f = open(os.path.join(path, file[:-3] + "_" + bfile[:-4] + "_peaks.bed"), "w")
    for i, region in enumerate(interval):
        print(region[0], region[1], region[2], end='\n', sep='\t', file=f)
    f.close()

    print("End")


@ray.remote
def filtering_tsv(filename, barcode_num, inte):
    temp_counts = np.zeros((len(barcode_num), len(inte)), dtype=int)
    sam_file = pysam.TabixFile(filename)

    for ii_, b_ in enumerate(inte):
        try:
            for read in sam_file.fetch(b_[0], b_[1], b_[2]):
                rows = read.split('\t')
                if barcode_num.__contains__(rows[3]):
                    if (int(rows[1]) > int(b_[1])) and (int(rows[1]) < int(b_[2])):
                        temp_counts[barcode_num[rows[3]], ii_] += 1
                    if (int(rows[2]) > int(b_[1])) and (int(rows[2]) < int(b_[2])):
                        temp_counts[barcode_num[rows[3]], ii_] += 1

        except:
            print("Failed to read: {}.".format(b_))

    return temp_counts


# COUNT FRAGMENTS IN PEAKS
def filtering_fragments(fragments, barcodes):
    print("Reading Barcodes")
    print(barcodes)
    path_barc, barc_file = os.path.split(barcodes)
    barc_d = []
    with open(barcodes, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            if (len(row)) > 0:
                barc_d.append(row[0])
    barc_d = np.array(barc_d)

    # FILTER TSV
    print("Filtering")
    print(fragments)
    with open(fragments) as f:
        with open(fragments + "_" + barc_file + '.tsv', 'w') as f_out:
            reader = csv.reader(f)
            for row in reader:
                if row[0].split('\t')[3] in barc_d:
                    f_out.write(row[0] + "\n")
