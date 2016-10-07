import unittest
import seqPipeline

class seqPipelineTest(unittest.TestCase):
	pipeline = seqPipeline()
	fastqFile = os.path.join('testFiles', 'fastqExample.fastq')
	samFile = os.path.join('testFiles', 'samExample.sam')
	def test_sam2bam2sortedBam2(samFile):
		

if __name__ == "__main__":
	unittest.main()
