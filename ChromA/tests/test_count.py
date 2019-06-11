from unittest import TestCase
from scipy.io import mmread
import numpy as np
import ChromA
import os


class TestChromA(TestCase):
    def test_count_correctly(self):
        # Package Directory
        dh = ChromA.data_handle
        pkg_dir = ChromA.__path__

        # Define Test Files and Directories
        tst_dt_dir = pkg_dir + "/test_data/count"
        files = [os.path.join(tst_dt_dir, "test_fragments.tsv"),
                 os.path.join(tst_dt_dir, "test_peaks.bed"), os.path.join(tst_dt_dir, "test_singlecell.csv")]
        order = [0, 1, 2]
        name = tmpdir.join("test_calc.log")

        # Create Logger and Usual Validations
        dh.build_logger(options.verbose_level, filename=name)
        logger = logging.getLogger()
        reads, inte = dh.count_validate_inputs(files=files, order=order)

        # Count cells in peaks
        logger.info("Counting Files: {}, {}".format(files[0], files[1]))
        dh.count_count(files=options.count, order=order, reads_array=reads, intervals=inte)

        # Load Result Files and Test or Correctness
        # import filecmp
        # filecmp.cmp('file1.bed', 'file1.bed')
        # filecmp.cmp('file1.tsv', 'file1.tsv')
        # filecmp.cmp('file1.mtx', 'file1.mtx')
