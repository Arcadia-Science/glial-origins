from string_functions import *

# Class definitions to handle biofiles between scripts
class sample_dict:
    def __init__(self, species, conditions, directory):
        self.species = species
        self.conditions = conditions
        self.directory = directory

# Class biofile is used to carry metadata for each file
class biofile:
    def __init__(self, filename, sample_dict):
        self.filename = filename
        self.species = sample_dict.species
        self.species_prefix = prefixify(self.species)
        self.directory = sample_dict.directory
        self.path = self.directory + self.filename
        self.filetype = self.filename.split('.')[-1]
        self.s3uri = None
    
    def push_to_s3(self):
        import subprocess
        subprocess(['aws', 's3', 'cp', self.directory + self.filename, self.s3uri])
        return None
    
    def get_from_s3(self, directory = ''):
        import os
        import subprocess
        
        if directory != '':
            self.directory = directory
        
        if not os.path.exists(directory + filename):
            subprocess(['aws', 's3', 'cp', self.s3uri, self.directory + self.filename])
        return None
    
    def add_s3uri(self, s3uri):
        if self.s3uri == None:
            self.s3uri = s3uri
        else:
            raise Exception('This file already has an S3 URI at ' + self.s3uri)
        return None

class genome_fasta_file(biofile):
    def __init__(self, filename, sample_dict):
        biofile.__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename

class transdecoder_cdna_file(biofile):
    def __init__(self, filename, sample_dict):
        biofile.__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/transcriptome/' + self.filename

class genome_gff_file(biofile):
    def __init__(self, filename, sample_dict):
        biofile.__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename

class genome_gtf_file(biofile):
    def __init__(self, filename, sample_dict):
        biofile.__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename

class gxc_file(biofile):
    def __init__(self, filename, sample_dict):
        biofile.__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename

        
        
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
    import subprocess
    filename = s3uri.split('/')[-1]
    
    directory = sample_dict.directory
    subprocess(['aws', 's3', 'cp', s3uri, directory + filename])
    
    outfile = fileclass(filename, sample_dict)
    outfile.add_s3uri(s3uri)
    
    return outfile

def gffread_gff_to_gtf(genome_gff_file, sample_dict, gffread_loc):
    import subprocess

    filename = genome_gff_file.filename.replace('.gff', '.gtf')
    output = genome_gtf_file(filename, sample_dict)
    
    subprocess.run([gffread_loc, genome_gff_file.path, '-T', '-o', output.path])
    
    return output