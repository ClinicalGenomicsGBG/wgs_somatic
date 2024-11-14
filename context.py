import os

class RunContext:
    def __init__(self, run_path):
        self.run_path = run_path
        self.run_name = os.path.basename(run_path)
        self.run_date = self.run_name.split('_')[0]
        self.run_flowcell = self.run_name.split('_')[-1].split('+')[0]
        self.run_tag = '_'.join([self.run_date, self.run_flowcell])
        #self.samplesheet_path = os.path.join(run_path, 'SampleSheet.csv')
        self.demultiplex_summary_path = os.path.join(run_path, 'demuxer.json')

        self.sample_contexts = []

    @property
    def demultiplex_complete(self):
        return os.path.exists(self.demultiplex_summary_path)

    def add_sample_context(self, Sctx):
        self.sample_contexts.append(Sctx)



class SampleContext:
    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.sample_name = self.sample_id.split('_')[0]
        self.slims_info = {}
        self.fastqs = []

    def add_fastq(self, paths):
        """Add fastq paths."""
        for path in paths:
            self.fastqs.append(path)
