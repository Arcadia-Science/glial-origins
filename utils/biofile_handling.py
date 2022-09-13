from string_functions import *

# Class definitions to handle biofiles between scripts
class sample_dict:
    def __init__(self, species: str, conditions: str, directory: str):
        self.species = species
        self.conditions = conditions
        self.directory = directory

# Class biofile is used to carry metadata for each file
class biofile:
    def __init__(self, filename: str, sample_dict: sample_dict, **kwargs):
        self.filename = filename
        self.species = sample_dict.species
        self.conditions = sample_dict.conditions
        self.species_prefix = prefixify(self.species)
        self.directory = sample_dict.directory
        self.path = self.directory + self.filename
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

class genome_fasta_file(biofile):
    def __init__(self, filename: str, sample_dict: sample_dict, version, **kwargs):
        biofile.__init__(self, filename, sample_dict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
        self.version = version

class transdecoder_cdna_file(biofile):
    def __init__(self, filename: str, sample_dict: sample_dict, genome_fasta_file, genome_annot_file, **kwargs):
        biofile.__init__(self, filename, sample_dict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/transcriptome/' + self.filename
        self.reference_genome = genome_fasta_file
        self.reference_annot = genome_annot_file

class genome_gff_file(biofile):
    def __init__(self, filename: str, sample_dict: sample_dict, genome_fasta_file, **kwargs):
        biofile.__init__(self, filename, sample_dict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = genome_fasta_file
    
    def to_gtf(self, GFFREAD_LOC: str):
        import subprocess
        
        self_sample = sample_dict(self.species, self.conditions, self.directory)

        filename = self.filename.replace('.gff', '.gtf')
        output = genome_gtf_file(filename, sample_dict = self_sample, genome_fasta_file = self.reference_genome)
        subprocess.run([GFFREAD_LOC, self.path, '-T', '-o', output.path])
        
        return output

class genome_gtf_file(biofile):
    def __init__(self, filename: str, sample_dict: sample_dict, genome_fasta_file, **kwargs):
        biofile.__init__(self, filename, sample_dict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = genome_fasta_file

class gxc_file(biofile):
    def __init__(self, filename: str, sample_dict: sample_dict, genome_fasta_file, genome_annot_file, **kwargs):
        biofile.__init__(self, filename, sample_dict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.reference_genome = genome_fasta_file
        self.reference_annot = genome_annot_file

class idmm_file(biofile):
    def __init__(self, filename: str, sample_dict: sample_dict, kind: str, sources: list, **kwargs):
        biofile.__init__(self, filename, sample_dict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        self.kind = kind
        
        
        
class docket:
    def __init__(self, sample_dict):
        self.species = sample_dict.species
        self.conditions = sample_dict.conditions
        self.directory = sample_dict.directory
        self.files = dict()
    
    def add_file(self, biofile):
        self.files[biofile.filename] = biofile
    
    def add_files(self, list):
        for biofile in list:
            self.files[biofile.filename] = biofile
    
    def add_keyfile(self, biofile, key):
        setattr(self, key, biofile)
        
    def add_keyfiles(self, dictionary):
        for key in dictionary:
            setattr(self, key, dictionary[key])
    
    def remove_file(self, key):
        del self.files[key]




# Functions for handling biofiles
        
# function for downloading and unzipping a file based on url
def url_download_biofile(url = str(), protocol = str(),
                         sample_dict = sample_dict, fileclass = biofile):
    
    raise Warning('deprecated! use biofile.get_from_url() method instead.')
    
    import os
    import subprocess
    
    species = sample_dict.species
    conditions = sample_dict.conditions
    directory = sample_dict.directory
    
    protocols = ['curl']
    if protocol not in protocols:
        raise ValueError("Invalid protocol. Expected one of: %s" % protocols)
        
    compression = url.split('.')[-1]
    
    compressions = ['gz', 'gzip', 'zip']
    if compression not in compressions:
        compression = 'none'
    
    filename = url.split('/')[-1]
    
    # make sure output folder ends with a backslash
    if directory[-1] != '/':
        outdir = directory + '/'
    else:
        outdir = directory
    output_loc = outdir + filename
    
    if compression != 'none':
        unzip_filename = filename.replace('.' + compression, '')
        unzip_loc = outdir + unzip_filename
        if os.path.exists(unzip_loc):
            return fileclass(unzip_filename, sample_dict) 
        elif os.path.exists(output_loc):
            _
        else:
            if protocol == 'curl':
                subprocess.run([protocol, url, '--output', output_loc])
                
        if compression == 'gz' or 'gzip':
            subprocess.run(['gunzip', output_loc])
            
        return fileclass(unzip_filename, sample_dict)
    else:
        if protocol == 'curl':
                subprocess.run([protocol, url, '--output', output_loc])
        return fileclass(filename, sample_dict)
        
def biofile_from_s3(s3uri, sample_dict, fileclass = biofile):
    
    raise Warning('deprecated! use biofile.get_from_s3() method instead.')
    
    import subprocess
    filename = s3uri.split('/')[-1]
    
    directory = sample_dict.directory
    subprocess(['aws', 's3', 'cp', s3uri, directory + filename])
    
    outfile = fileclass(filename, sample_dict)
    outfile.add_s3uri(s3uri)
    
    return outfile

def gffread_gff_to_gtf(genome_gff_file, sample_dict, gffread_loc):
    
    raise Warning('deprecated! use genome_gff_file.to_gtf() method instead.')
    
    import subprocess

    filename = genome_gff_file.filename.replace('.gff', '.gtf')
    output = genome_gtf_file(filename, sample_dict)
    
    subprocess.run([gffread_loc, genome_gff_file.path, '-T', '-o', output.path])
    
    return output