from string_functions import *

# Class definitions to handle BioFiles between scripts
class SampleDict:
    def __init__(self, species: str, conditions: str, directory: str):
        self.species = species
        self.conditions = conditions
        self.directory = directory

class Docket:
    def __init__(self, SampleDict):
        self.species = SampleDict.species
        self.conditions = SampleDict.conditions
        self.directory = SampleDict.directory
        self.files = dict()
        
    def add_file(self, BioFile):
        self.files[BioFile.filename] = BioFile
        
    def add_files(self, list):
        for BioFile in list:
            self.files[BioFile.filename] = BioFile
            
    def add_keyfile(self, BioFile, key):
        setattr(self, key, BioFile)
        
    def add_keyfiles(self, dictionary):
        for key in dictionary:
            setattr(self, key, dictionary[key])
            
    def remove_file(self, key):
        del self.files[key]
        
        
# Class BioFile is used to carry metadata for each file
class BioFile:
    def __init__(self, filename: str, SampleDict: SampleDict, **kwargs):
        self.filename = filename
        self.species = SampleDict.species
        self.conditions = SampleDict.conditions
        self.directory = SampleDict.directory
        self.path = self.directory + self.filename
        self.species_prefix = prefixify(self.species)
        self.s3uri = None
        
        if 'url' and 'protocol' in kwargs.keys():
            self.get_from_url(kwargs['url'], kwargs['protocol'])
        elif 's3uri' in kwargs.keys():
            self.get_from_s3(kwargs['s3uri'])
    
    def add_s3uri(self, s3uri: str):
        if self.s3uri == None:
            self.s3uri = s3uri
        else:
            raise Exception('This file already has an S3 URI at ' + self.s3uri)
        return None
    
    def get_from_s3(self, s3uri: str, directory = ''):
        import os
        import subprocess
        
        if self.s3uri == None:
            self.s3uri = s3uri
            self.filetype = self.s3uri.split('.')[-1]
        else:
            raise Exception('This file already has an S3 URI at ' + self.s3uri)
        
        if directory != '':
            self.directory = directory
        
        if not os.path.exists(directory + filename):
            subprocess(['aws', 's3', 'cp', self.s3uri, self.directory + self.filename])
        return self
    
    def push_to_s3(self):
        import subprocess
        
        if self.s3uri is None:
            pass
        
        subprocess(['aws', 's3', 'cp', self.directory + self.filename, self.s3uri])
        return None

    def get_from_url(self, url: str, protocol: str):
        import os
        import subprocess
        
        self.url = url

        protocols = ['curl']
        if protocol not in protocols:
            raise ValueError("Invalid protocol. Expected one of: %s" % protocols)   
        
        compression = self.url.split('.')[-1]
        compressions = ['gz', 'gzip', 'zip']
        
        if compression not in compressions:
            compression = 'none'
        
        filename = self.url.split('/')[-1]
        
        if compression != 'none':
            output_loc = self.directory + filename
            self.filename = filename.replace('.' + compression, '')
            self.path = self.directory + self.filename
            self.filetype = self.filename.split('.')[-1]
            
            if os.path.exists(self.path):
                return self
            elif os.path.exists(output_loc):
                _
            else:
                if protocol == 'curl':
                    subprocess.run([protocol, url, '--output', output_loc])
                
                if compression == 'gz' or 'gzip':
                    subprocess.run(['gunzip', output_loc])
                return self
        else:
            self.filename = filename
            self.path = self.directory + self.filename
            self.filetype = self.filename.split('.')[-1]
            
            if protocol == 'curl':
                subprocess.run([protocol, url, '--output', self.path])
            return self

class GenomeFastaFile(BioFile):
    def __init__(self, filename: str, SampleDict: SampleDict, version: str, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.version = version
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
        
    def get_transdecoder_cdna_gtf(self, GenomeGtfFile, TRANSDECODER_LOC, **kwargs):
        PERL_SCRIPT_LOC = TRANSDECODER_LOC + 'util/gtf_genome_to_cdna_fasta.pl'
        
        import os
        import subprocess
        
        self_SampleDict = SampleDict(self.species, self.conditions, self.directory)
        
        cdna_name = self.filename.replace('.' + self.filetype, '_cDNA.' + self.filetype)
        cdna_output = TransdecoderCdnaFile(
            filename = cdna_name,
            SampleDict = self_SampleDict,
            GenomeFastaFile = self,
            GenomeAnnotFile = GenomeGtfFile
        )
        if os.path.exists(cdna_output.path):
            return cdna_output
        else:
            output_file = open(cdna_output.path, "w")
            
            subprocess.call([PERL_SCRIPT_LOC, GenomeGtfFile.path, self.path], stdout=output_file)
            return cdna_output


class TransdecoderCdnaFile(BioFile):
    def __init__(self, filename: str, SampleDict: SampleDict, GenomeFastaFile, GenomeAnnotFile, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/transcriptome/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile

    def to_pep_files(self, TDLONGORF_LOC, TDPREDICT_LOC):
        import subprocess
        
        suffixes = ['.bed', '.cds', '.gff3', '.pep']
        temp_dir = self.directory + ('_'.join([self.species, self.conditions])) + '.transdecoder_dir/'
        
        subprocess.run([TDLONGORF_LOC, '-t', self.path, '-O', temp_dir])
        subprocess.run([TDPREDICT_LOC, '-t', self.path, '-O', temp_dir])
        
        output_dict = {}
        self_SampleDict = SampleDict(self.species, self.conditions, self.directory)
        
        for suffix in suffixes:
            output_filename = self.filename + '.transdecoder' + suffix
            output_dict['transdecoder_' + suffix.replace('.', '')] = TransdecoderOutFile(
                output_filename, 
                self_SampleDict, 
                GenomeFastaFile = self.reference_genome, 
                GenomeAnnotFile = self.reference_annot,
                TransdecoderCdnaFile = self
            )
        
        return output_dict
    
class TransdecoderOutFile(BioFile):
    def __init__(self, filename, SampleDict, GenomeFastaFile, GenomeAnnotFile, TransdecoderCdnaFile, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/proteome/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile
        self.reference_cDNA = TransdecoderCdnaFile

class GenomeGffFile(BioFile):
    def __init__(self, filename: str, SampleDict: SampleDict, GenomeFastaFile: GenomeFastaFile, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = GenomeFastaFile
    
    def to_gtf(self, GFFREAD_LOC):
        import subprocess
        
        self_SampleDict = SampleDict(self.species, self.conditions, self.directory)

        filename = self.filename.replace('.gff', '.gtf')
        output = GenomeGtfFile(filename, SampleDict = self_SampleDict, GenomeFastaFile = self.reference_genome)
        subprocess.run([GFFREAD_LOC, self.path, '-T', '-o', output.path])
        
        return output

class GenomeGtfFile(BioFile):
    def __init__(self, filename, SampleDict, GenomeFastaFile, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = GenomeFastaFile

class GxcFile(BioFile):
    def __init__(self, filename: str, SampleDict: SampleDict, GenomeFastaFile, GenomeAnnotFile, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile

class IdmmFile(BioFile):
    def __init__(self, filename, SampleDict, kind, sources, **kwargs):
        super().__init__(filename = filename, SampleDict = SampleDict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        self.kind = kind