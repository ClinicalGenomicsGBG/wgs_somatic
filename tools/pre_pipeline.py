from collections import defaultdict

from slims import Sample

def pre_pipeline(runs: defaultdict[str, list[Sample]]):
    """ 
    pre-pipeline which handles:
    - check presence of files locally and if not, download them from the server
    - merge fastq for normals and tumors
    - basic stats and expected coverage
    """
    
