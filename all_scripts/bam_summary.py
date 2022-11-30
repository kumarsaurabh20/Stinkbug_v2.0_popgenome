import os
import pysam


class Samstats:

    def __init__(self, bam):
        self.bam = bam

    def sam_check(self):
        # name_tuple = os.path.splitext(os.path.basename(self.bam))
        if os.path.exists(self.bam):
            if os.path.exists(self.bam + 'bai'):
                return True
            else:
                pysam.index(self.bam)
                return True
        else:
            print("BAM file was not found!")
            return False

    def mapped_reads(self, position=""):  # position = 'ref 10 100'
        bamfile = pysam.AlignmentFile(self.bam, 'rb')
        try:
            bamfile.header['HD']['SO']
        except KeyError:
            print("BAM file is not sorted!!")
            exit()
            # pysam.sort('-@', 8, '-o', name_tuple[0] + 'sorted.bam', self.bam)
            # os.remove(self.bam)
            # os.rename(name_tuple[0] + 'sorted.bam', self.bam)
        if not position:
            count = bamfile.count()
        else:
            count = bamfile.count(position)
        return count
        bamfile.close()
