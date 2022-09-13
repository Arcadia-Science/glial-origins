import os
import subprocess

from string_functions import *

# Class definitions to handle biofiles between scripts
class sample_dict:
    def __init__(self, species, conditions, directory):
        self.species = species
        self.conditions = conditions
        self.directory = directory

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

# Class BioFile is used to carry metadata for each file
class BioFile:
    def __init__(self, filename, sample_dict):
        self.filename = filename
        self.species = sample_dict.species
        self.directory = sample_dict.directory
        self.s3uri = None

    @property
    def path(self):
        return self.directory + self.filename

    @property
    def filetype(self):
        return self.filename.split('.')[-1]

    @property
    def species_prefix(self):
        return prefixify(self.species)

    def push_to_s3(self):
        if self.s3uri is None:
            pass
        subprocess(['aws', 's3', 'cp', self.path, self.s3uri])

class GenomeFastaFile(BioFile):
    def __init__(self, filename, sample_dict, version):
        super().__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
        self.version = version

class TransdecoderCdnaFile(BioFile):
    def __init__(self, filename, sample_dict, genome_fasta_file, genome_annot_file):
        super().__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
        self.reference = genome_fasta_file.filename
        self.reference_version = genome_fasta_file.version
        self.annotation = genome_annot_file.filename

class GenomeGffFile(BioFile):
    def __init__(self, filename, sample_dict, genome_fasta_file):
        super().__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference = genome_fasta_file.filename
        self.version = genome_fasta_file.version

class GenomeGtfFile(BioFile):
    def __init__(self, filename, sample_dict, genome_fasta_file):
        super().__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference = genome_fasta_file.filename
        self.version = genome_fasta_file.version

class GxcFile(BioFile):
    def __init__(self, filename, sample_dict, genome_fasta_file):
        super().__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.reference = genome_fasta_file.filename
        self.reference_version = genome_fasta_file.version

class IdmmFile(BioFile):
    def __init__(self, filename, sample_dict, fields):
        super().__init__(self, filename, sample_dict)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.fields = fields

# Functions for handling biofiles

# function for downloading and unzipping a file based on url
def url_download_biofile(url = str(), protocol = str(),
                         sample_dict = sample_dict, fileclass = biofile, **kwargs):
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
            return fileclass(unzip_filename, sample_dict, **kwargs)
        elif os.path.exists(output_loc):
            _
        else:
            if protocol == 'curl':
                subprocess.run([protocol, url, '--output', output_loc])

        if compression == 'gz' or 'gzip':
            subprocess.run(['gunzip', output_loc])

        return fileclass(unzip_filename, sample_dict, **kwargs)
    else:
        if protocol == 'curl':
                subprocess.run([protocol, url, '--output', output_loc])
        return fileclass(filename, sample_dict, **kwargs)

def gffread_gff_to_gtf(genome_gff_file, sample_dict, gffread_loc):
    filename = genome_gff_file.filename.replace('.gff', '.gtf')
    output = genome_gtf_file(filename, sample_dict)

    subprocess.run([gffread_loc, genome_gff_file.path, '-T', '-o', output.path])

    return output
