from string_functions import *
from install_locs import *
import os, subprocess, re
from typing import Union

if not os.path.exists(GLOBAL_OUTPUT_DIRECTORY):
    os.mkdir(GLOBAL_OUTPUT_DIRECTORY)
    print(GLOBAL_OUTPUT_DIRECTORY, 'does not exist; making it now')

def s3_transfer(to_loc: str, from_loc: str):
    """Transfers files to and from AWS S3 and local.  
    
    One of the two parameters must be an AWS S3 URI in string format.  
    
    Args:
        to_loc (str): path of the destination file.
        from_loc (str): path of the origin file.
    """
    subprocess.run(['aws', 's3', 'cp', to_loc, from_loc])

# Makes a directory based on species and conditions
# If species
def make_output_directory(species: str, conditions: str, stringonly = False):
    """Creates an output directory based on species and condition parameters.  
    
    Does not create a directory if it already exists.  
    Output directory is placed in the `GLOBAL_OUTPUT_DIRECTORY`, as defined by `env/install_locs.py`.
    Default `GLOBAL_OUTPUT_DIRECTORY` is `output/`.
    
    Args:
        species (str): species name in `Genus_species` format.
        conditions (str): unique conditions identifier for dataset.
        stringonly (bool, optional): whether to only output the string without creating a directory
    
    Returns:
        output_directory (str): path of the output directory, in format `Gspe_conditions`.
    """
    from string_functions import prefixify
    
    species_prefix = prefixify(species)

    # Specify folder as destination for file downloads
    output_directory = GLOBAL_OUTPUT_DIRECTORY + prefixify(species) + '_' + conditions + '/'
    
    if stringonly:
        return output_directory

    if not os.path.exists(output_directory):
        print('creating', output_directory)
        os.mkdir(output_directory)
    else:
        print(output_directory, 'already exists')
    
    return output_directory

class dummy_object(object):
    """simple dummy object to enable dot access to keys"""
    pass

# Class definitions to handle BioFiles between scripts
class SampleDict:
    """A dictionary containing dataset-specific fields: species, conditions, and directory.
    
    This class is used to uniquely associate BioFile objects with specific datasets.
    
    Args:
        species (str): species name in 'Genus_species' format.
        conditions (str): unique conditions identitier for dataset.
        directory (str): output directory for files in the dataset.
    """
    def __init__(self, species: str, conditions: str, directory: str):
        self.species = species
        self.conditions = conditions
        self.directory = directory

class BioFileDocket:
    """BioFileDocket objects collect BioFiles and relate them to one another.
    
    Important files called `keyfiles` get their own uniquely-named attribute.  
    These files can be accessed using a dot operator.
    For example, `BFD.genome_fasta` would return a GenomeFastaFile object.
    `BFD.genome_fasta.path` would return the path to that GenomeFastaFile object.
    
    Args:
        species (str): species name in 'Genus_species' format.
        conditions (str): unique conditions identitier for dataset.
    
    Attributes:
        directory (str): output directory for files in the dataset.
    """
    def __init__(self, species: int, conditions: int):
        if '_' not in species:
            raise Exception('Please include an underscore in the species name with format:\nGenus_species')
        self.species = species
        
        if any(not c.isalnum() for c in conditions):
            raise Exception('Conditions can only include alphanumeric characters')
        self.conditions = conditions
        
        directory = make_output_directory(self.species, self.conditions)
        self.directory = directory
        print('Files will be saved into', self.directory)
        
        self.files = dict()
        self.metadata = dummy_object()
    
    @property
    def s3uri(self):
        """str: the S3 URI of the BioFileDocket."""
        return 's3://arcadia-reference-datasets/glial-origins-pkl/' + self.dill_filename
    
    @property
    def dill_filename(self):
        """str: the unique filename of the BioFileDocket .pkl file."""
        return '_'.join([prefixify(self.species), self.conditions, 'sample_BioFileDocket.pkl'])
    
    @property
    def dill_filepath(self):
        """str: the path of the BioFileDocket .pkl file."""
        return self.directory + self.dill_filename
    
    @property
    def sampledict(self):
        """:obj:`SampleDict` the full SampleDict of the BioFileDocket."""
        return SampleDict(self.species, self.conditions, self.directory)
    
    def set_metadata(self, key: str, value: any, replace = False):
        """Adds a metadata feature using a unique key identifier.
        
        Checks to make sure that key does not already exist.
        Also makes sure key is only alphanumeric or underscores.
        
        Args:
            key (str): a unique key
            value (obj): any field you want to record
        """
        if not re.match(r'^\w+$', key):
            raise Exception('key can only include alphanumeric characters and underscores')
        if hasattr(self.metadata, key):
            if replace:
                print('overwriting ' + key)
            else:
                print('key "' + key + '" already exists, ignoring')
                return
        setattr(self.metadata, key, value)
       
    def set_taxid(self, taxid: Union[str, int]):
        """Adds a taxid attribute to the BioFileDocket.
        
        Tolerates either `int` or `str` input.
        """
        self.metadata.taxid = str(taxid)
        
    def add_file(self, BioFile):
        """Places a BioFile object into the `files` attribute.
        
        The objects can be accessed using their filename.
        
        Note:
            The actual files of BioFile objects in this list are not automatically uploaded.
        """
        self.files[BioFile.filename] = BioFile
        
    def add_files(self, list):
        """Places a list of files into the `files` attribute."""
        for BioFile in list:
            self.add_file(BioFile)
    
    def remove_file(self, filename):
        """Removes a file based on filename from the files attribute."""
        del self.files[filename]
            
    def add_keyfile(self, key: str, BioFile):
        """Adds a BioFile object using a unique key identifier.
        
        Checks to make sure that key does not already exist.
        Also makes sure key is only alphanumeric or underscores.
        
        Args:
            key (str): a unique key
            BioFile (:obj:`BioFile`): a BioFile object
        """
        if not re.match(r'^\w+$', key):
            raise Exception('key can only include alphanumeric characters and underscores')
        if hasattr(self, key):
            raise Exception('key "' + key + '" already exists')
        setattr(self, key, BioFile)
        
    def add_keyfiles(self, dictionary: dict):
        """Adds a dictionary of BioFile objects using key:BioFile object pairs."""
        for key in dictionary:
            self.add_keyfile(key, dictionary[key])
    
    def remove_keyfile(self, key, warn = True):
        """Deletes a keyfile from the BioFileDocket if warn == False."""
        if not warn:
            delattr(self, key)
        else:
            raise Warning('If you want to delete this keyfile, set warn = False.')
    
    def pickle(self):
        """Creates a .pkl file for the BioFileDocket object.
        
        The filename and path are automatically generated.
        """
        import dill
        
        with open(self.dill_filepath, 'wb') as file:
            dill.dump(self, file)
    
    def unpickle(self):
        """Unpickles a .pkl file for the BioFileDocket object.
        
        The filename and path are automatically generated.
        
        Returns:
            self :obj:`BioFileDocket`: returns the unpickled BioFileDocket object.
        
        Raises:
            FileNotFoundError: if the .pkl file doesn't exist yet.
        """
        
        import dill
        import os
        
        if not os.path.exists(self.dill_filepath):
            raise FileNotFoundError("Can't unpickle file; .pkl file doesn't exist yet at " + self.dill_filepath)

        with open(self.dill_filepath, 'rb') as file:
            self = dill.load(file)
        
        return self
    
    def local_to_s3(self, overwrite = False):
        """Uploads all keyfiles that are BioFiles to S3 using their s3uri attributes.
        
        Args:
            overwrite (bool): decide whether to overwrite existing files. Defaults to False.
        """
        files = {i:j for i,j in dict(vars(self)).items() if isinstance(j, BioFile)}

        for file in files.values():
            file.push_to_s3(overwrite)
    
    def s3_to_local(self, overwrite = False):
        """Downloads all keyfiles that are BioFiles from S3 using their s3uri attributes.
        
        Args:
            overwrite (bool): decide whether to overwrite existing files. Defaults to False.
        """
        files = {i:j for i,j in dict(vars(self)).items() if isinstance(j, BioFile)}
        
        for file in files.values():
            file.get_from_s3(overwrite)
    
    def get_from_s3(self, overwrite = False):
        """Downloads the .pkl file for the BioFileDocket from AWS S3.
        
        Args:
            overwrite (bool): decide whether to overwrite existing files. Defaults to False.
        """
        import os
        import subprocess
        
        if not os.path.exists(self.dill_filepath):
            s3_transfer(self.s3uri, self.dill_filepath)
            return self
        elif overwrite:
            s3_transfer(self.s3uri, self.dill_filepath)
            return self
        else:
            print('file', self.dill_filename, 'already exists at', self.dill_filepath)
            return self
    
    def push_to_s3(self, overwrite = False):
        """Uploads the .pkl file for the BioFileDocket to AWS S3.
        
        Args:
            overwrite (bool): decide whether to overwrite existing files. Defaults to False.
        """
        import subprocess
        
        bucket = self.s3uri.lstrip('s3://').split('/')[0]
        file_key = self.s3uri.split(bucket)[1].lstrip('/')
        
        # check if file exists in S3 already
        output = subprocess.run(['aws', 's3api', 'head-object', '--bucket', bucket, '--key', file_key], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        if 'Not Found' in str(output.stderr):
            s3_transfer(self.dill_filepath, self.s3uri)
        elif overwrite:
            print(self.dill_filename, 'already exists in S3 bucket; overwriting.')
            s3_transfer(self.dill_filepath, self.s3uri)
        else:
            print(self.dill_filename, 'already exists in S3 bucket, skipping upload. set overwrite = True to overwrite the existing file.')
        return None

class MultiSpeciesBioFileDocket:
    def __init__(self, species_dict: dict, global_conditions: str, analysis_type):
        self.species_dict = species_dict
        self.global_conditions = global_conditions
        self.analysis_type = analysis_type
        self.SampleDicts = {}
        self.metadata = dummy_object()
        
        # create dummy species_SampleDicts
        for species in species_dict:
            conditions = species_dict[species]
            species_prefix = prefixify(species)
            output_folder = make_output_directory(species, conditions, stringonly = True)
            species_SampleDict = SampleDict(species, conditions, output_folder)
            self.SampleDicts[species_prefix] = species_SampleDict
        
        directory = make_output_directory(self.species_concat, '_'.join([self.global_conditions, self.analysis_type]))
        self.directory = directory
    
    @property
    def species_concat(self):
        return ''.join(sorted(self.SampleDicts.keys()))
    
    def get_BioFileDockets(self):
        import dill
        
        self.species_BioFileDockets = {}
        
        for sampledict in self.SampleDicts.values():
            species_prefix = prefixify(sampledict.species)
            
            dill_filename = '_'.join([species_prefix, sampledict.conditions, 'sample_BioFileDocket.pkl'])
            dill_filepath = sampledict.directory + dill_filename
            dill_s3uri = 's3://arcadia-reference-datasets/glial-origins-pkl/' + dill_filename
            
            if not os.path.exists(dill_filepath):
                s3_transfer(dill_s3uri, dill_filepath)
            
            with open(dill_filepath, 'rb') as file:
                self.species_BioFileDockets[species_prefix] = dill.load(file)
        return None
    
    def s3_to_local(self, overwrite = False):
        for species_prefix in self.species_BioFileDockets:
            self.species_BioFileDockets[species_prefix].s3_to_local(overwrite)
    
    def local_to_s3(self, overwrite = False):
        for species_prefix in self.species_BioFileDockets:
            self.species_BioFileDockets[species_prefix].local_to_s3(overwrite)
        
        
# Class BioFile is used to carry metadata for each file
class BioFile:
    def __init__(self, sampledict: SampleDict, filename = '', url = None, protocol = None, unzip = True, s3uri = None):
        self.filename = filename
        self.species = sampledict.species
        self.conditions = sampledict.conditions
        self.directory = sampledict.directory
        self.s3uri = s3uri
        
        if url != None and protocol != None:
            self.get_from_url(url = url, protocol = protocol, filename = filename, unzip = unzip)
        elif s3uri != None:
            self.get_from_s3(s3uri)
        elif filename == '':
            raise Exception('File must minimally have a filename.\nTo download from url, pass both a url and protocol.\nTo download from s3uri, pass s3uri.')
            
    @property
    def exists(self):
        return os.path.exists(self.path)
    
    @property
    def path(self):
        return self.directory + self.filename
    
    @property
    def species_prefix(self):
        if '_' not in self.species:
            return self.species
        else:
            return prefixify(self.species)
    
    @property
    def sampledict(self):
        return SampleDict(self.species, self.conditions, self.directory)
    
    @property
    def filetype(self):
        return self.filename.split('.')[-1]
    
    def add_s3uri(self, s3uri: str):
        if self.s3uri == None:
            self.s3uri = s3uri
        else:
            raise Exception('This file already has an S3 URI at ' + self.s3uri)
        return None
    
    def get_from_s3(self, overwrite = False):
        import os
        import subprocess
        
        if self.filename == '':
            self.filename = self.s3uri.split('/')[-1]
        
        if not os.path.exists(self.path):
            s3_transfer(self.s3uri, self.path)
            return self
        elif overwrite:
            s3_transfer(self.s3uri, self.path)
            return self
        else:
            print('file', self.filename, 'already exists at', self.path)
            return self
    
    def push_to_s3(self, overwrite = False):
        import subprocess
        
        if self.s3uri is None:
            raise Exception(self.path, 'has no s3uri. Use add_s3uri() to add a URI.')
        
        bucket = self.s3uri.lstrip('s3://').split('/')[0]
        file_key = self.s3uri.split(bucket)[1].lstrip('/')
        
        # check if file exists in S3 already
        # suppress stdout and save stderr as file to detect non-existing files
        output = subprocess.run(['aws', 's3api', 'head-object', '--bucket', bucket, '--key', file_key], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        if 'Not Found' in str(output.stderr):
            s3_transfer(self.path, self.s3uri)
        elif overwrite:
            print(self.filename, 'already exists in S3 bucket; overwriting.')
            s3_transfer(self.path, self.s3uri)
        else:
            print(self.filename, 'already exists in S3 bucket, skipping upload. set overwrite = True to overwrite the existing file.')
        
        return None

    def get_from_url(self, url: str, protocol: str, filename = '', unzip = True):
        import os
        import subprocess
        
        self.url = url

        protocols = ['curl', 'wget', 'rsync']
        if protocol not in protocols:
            raise ValueError("Invalid protocol. Expected one of: %s" % protocols)   
        
        compression = self.url.split('.')[-1]
        compressions = ['gz', 'gzip', 'zip']
        
        if compression not in compressions:
            compression = 'none'
            
        if filename == '':
            filename = self.url.split('/')[-1]
        
        if compression != 'none':
            output_loc = self.directory + filename
            self.filename = filename.replace('.' + compression, '')
            
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
            
            if os.path.exists(self.path):
                return self
            
            if protocol == 'curl':
                subprocess.run([protocol, url, '--output', self.path])
            return self

class GenomeFastaFile(BioFile):
    def __init__(self, sampledict: SampleDict, version: str, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.version = version
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
    
    def rename_RefSeq_chromosomes(self, replace = False):
        from Bio import SeqIO

        new_filename = self.filename + '.renamed.fa'
        new_filepath = self.directory + new_filename

        with open(self.path) as original, open(new_filepath, 'w') as corrected:
            records = SeqIO.parse(self.path, 'fasta')
            for record in records:      
                if 'chromosome' in record.description:
                    print('starting name:', record.description)
                    new_id = str(record.description).split('chromosome')[1].split(',')[0].strip(' ')
                    record.id = new_id
                    record.description = new_id
                    print('new name:', record.id)
                SeqIO.write(record, corrected, 'fasta')
        
        if replace:
            self.filename = new_filename
            
            return None
        
        else:
            new_file = GenomeFastaFile(
                filename = new_filename,
                sampledict = self.sampledict,
                version = self.version
            )
            return new_file
        
    def get_transdecoder_cdna_gtf(self, GenomeGtfFile, TRANSDECODER_LOC, **kwargs):
        PERL_SCRIPT_LOC = TRANSDECODER_LOC + 'util/gtf_genome_to_cdna_fasta.pl'
        
        import os
        import subprocess
        
        cdna_name = self.filename.replace('.' + self.filetype, '_cDNA.' + self.filetype)
        cdna_output = TransdecoderCdnaFile(
            filename = cdna_name,
            sampledict = self.sampledict,
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
    def __init__(self, sampledict: SampleDict, GenomeFastaFile, GenomeAnnotFile, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
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
        
        for suffix in suffixes:
            output_filename = self.filename + '.transdecoder' + suffix
            output_file = TransdecoderOutFile(
                filename = output_filename, 
                sampledict = self.sampledict, 
                GenomeFastaFile = self.reference_genome, 
                GenomeAnnotFile = self.reference_annot,
                TransdecoderCdnaFile = self
            )
            subprocess.run(['mv', output_file.filename, output_file.path])
            output_dict['transdecoder_' + suffix.replace('.', '')] = output_file
        
        return output_dict
    
class TransdecoderOutFile(BioFile):
    def __init__(self, sampledict, GenomeFastaFile, GenomeAnnotFile, TransdecoderCdnaFile, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/proteome/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile
        self.reference_cDNA = TransdecoderCdnaFile

class GenomeGffFile(BioFile):
    def __init__(self, sampledict: SampleDict, GenomeFastaFile: GenomeFastaFile, filename = '', **kwargs):
        import subprocess
        import os
        
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        
        if self.filename.split('.')[-1] == 'gff3':
            new_filename = self.filename.replace('gff3', 'gff')
            new_path = self.directory + new_filename
            subprocess.run(['mv', self.path, new_path])
            self.filename = new_filename
        
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = GenomeFastaFile
    
    def to_gtf(self, GFFREAD_LOC, keep_all = False):
        import subprocess
        import os
        
        filename = self.filename.replace('.gff', '.gtf')
        output = GenomeGtfFile(filename = filename, sampledict = self.sampledict, GenomeFastaFile = self.reference_genome)
        
        if os.path.exists(output.path):
            print('Converted file', output.filename, 'already exists at:\n', output.path)
            return output
        
        subprocess.run([GFFREAD_LOC, self.path, '-T', '-o', output.path])
        
        return output

class GenomeGtfFile(BioFile):
    def __init__(self, sampledict, GenomeFastaFile, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = GenomeFastaFile

class CellRangerFileGroup:
    def __init__(self, sampledict: SampleDict, barcodes_url, features_url, matrix_url, protocol):
        self.sampledict = sampledict
        self.barcodes = CellRangerBarcodesFile(sampledict = sampledict, url = barcodes_url, protocol = protocol)
        self.matrix = CellRangerMatrixFile(sampledict = sampledict, url = matrix_url, protocol = protocol)
        self.features = CellRangerFeaturesFile(sampledict = sampledict, url = features_url, protocol = protocol)
    
    # Converts a group of CellRanger files into a single dense csv matrix
    def to_gxc(self, filename, GenomeFastaFile, GenomeAnnotFile, overwrite = False):
        import csv, gzip, os
        import scipy.io
        import pandas as pd
        
        # read in MEX format matrix as table
        mat = scipy.io.mmread(self.matrix.path)
         
        # list of transcript ids, e.g. 'ENSG00000243485'
        feature_ids = [row[0] for row in csv.reader(open(self.features.path), delimiter="\t")]
 
        # list of gene names, e.g. 'MIR1302-2HG'
        gene_names = [row[1] for row in csv.reader(open(self.features.path), delimiter="\t")]
 
        barcodes = [row[0] for row in csv.reader(open(self.barcodes.path), delimiter="\t")]
        
        # transform table to pandas dataframe and label rows and columns
        matrix = pd.DataFrame.sparse.from_spmatrix(mat)
        matrix.columns = barcodes
        matrix.insert(loc=0, column="gene_name", value=gene_names)

        # display matrix
        display(matrix.head(10))
        
        # save the table as a CSV (note the CSV will be a very large file)
        output = GxcFile(self.sampledict, GenomeFastaFile, GenomeAnnotFile, filename)
        
        if not os.path.exists(output.path):
            matrix.to_csv(output.path, index=False, sep = '\t')
        else:
            print('File already exists at', output.path, '. Set overwrite = True to overwrite.')
        
        return output

class CellRangerBarcodesFile(BioFile):
    def __init__(self, sampledict: SampleDict, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        
class CellRangerMatrixFile(BioFile):
    def __init__(self, sampledict: SampleDict, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename

class CellRangerFeaturesFile(BioFile):
    def __init__(self, sampledict: SampleDict, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename

class LoomFile(BioFile):
    def __init__(self, sampledict: SampleDict, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
    
    def to_gxc(self, filename, GenomeFastaFile, GenomeAnnotFile, overwrite = False):
        print('hi')
        
class GxcFile(BioFile):
    def __init__(self, sampledict: SampleDict, GenomeFastaFile, GenomeAnnotFile, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile

class ExcFile(BioFile):
    def __init__(self, sampledict: SampleDict, gxcfile, embedding, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.original = gxcfile
        self.embedding = embedding

class IdmmFile(BioFile):
    def __init__(self, sampledict, kind, sources: list, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        self.kind = kind
        
class CellAnnotFile(BioFile):
    def __init__(self, sampledict, sources: list, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.sources = sources

class UniprotIDMapperFile(IdmmFile):
    def __init__(self, sampledict, kind, sources: list, from_type, to_type, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, kind = kind, sources = sources, **kwargs)
        self.from_type = from_type
        self.to_type = to_type

class UniProtTaxidListFile(BioFile):
    def __init__(self, sampledict, taxid, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        
        self.taxid = str(taxid)
        
        if filename == '':
            self.filename = '_'.join([self.species_prefix, 'uniprot-taxid', self.taxid, 'genes.tsv'])
        
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        
        if not os.path.exists(self.path):
            with open(self.path, "w") as outfile:
                subprocess.run(['curl', 
                  'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_id%2Corganism_name&format=tsv&query=%28%28taxonomy_id%3A' + self.taxid + '%29%29'],
                  stdout = outfile)
        
class GeneListFile(BioFile):
    def __init__(self, sampledict, sources: list, genes: list, identifier: str, filename = '', **kwargs):
        from string_functions import make_gene_list, prefixify
         
        if filename == '':
            print('filename is ignored and generated by input as ', new_filename)
        else:
            new_filename = filename

        super().__init__(filename = new_filename, sampledict = sampledict, **kwargs)

        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        
        make_gene_list(genes, self.path)
    
    def get_uniprot_ids(self, ID_MAPPER_LOC, from_type, to_type):
        import subprocess
        
        num_lines = sum(1 for line in open(self.path))
        
        if num_lines > 100000:
            subprocess.run(['split', '-l', '100000', '-d', '--additional-suffix=_split_ids.txt', 
                            self.path, self.path.replace('_ids.txt', '')])
            for file in [i for i in os.listdir(self.directory) if 'split_ids' in i]:
                subprocess.run([ID_MAPPER_LOC, file, from_type, to_type])
            
        else:
            filename = self.filename.replace('_ids.txt', '_UniProtIDs.txt')
        
            subprocess.run([ID_MAPPER_LOC, self.path, from_type, to_type])
            output = UniprotIDMapperFile(
                filename = filename, sampledict = self.sampledict, kind = 'UniprotIDMapper', 
                sources = self.sources, from_type = from_type, to_type = to_type)
        
            return output
    
    def get_uniprot_ids_by_genename(self, ID_MAPPER_LOC, from_type, to_type, taxid):
        import subprocess
        
        filename = self.filename.replace('_ids.txt', '_UniProtIDs.txt')
        
        subprocess.run([ID_MAPPER_LOC, self.path, from_type, to_type, taxid])
        output = UniprotIDMapperFile(
            filename = filename, sampledict = self.sampledict, kind = 'UniprotIDMapper', 
            sources = self.sources, from_type = from_type, to_type = to_type)
        
        return output

class AlphaFoldFiles(BioFile):
    def __init__(self, sampledict: SampleDict, gxcfile, embedding, filename = '', **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.original = gxcfile
        self.embedding = embedding

class MultiSpeciesFile(BioFile):
    def __init__(self, multispeciesbiofiledocket, filename = ''):
        self.filename = filename
        self.multispeciesbiofiledocket = multispeciesbiofiledocket
        self.species = list(self.multispeciesbiofiledocket.species_dict.keys())
        self.global_conditions = self.multispeciesbiofiledocket.global_conditions
        self.directory = self.multispeciesbiofiledocket.directory
        self.species_concat = self.multispeciesbiofiledocket.species_concat
        self.s3uri = None
    
    @property
    def path(self):
        return self.directory + self.filename

class OrthoFinderOutputFile(MultiSpeciesFile):
    def __init__(self, multispeciesbiofiledocket, directory: str, filename = ''):
        super().__init__(filename = filename, multispeciesbiofiledocket = multispeciesbiofiledocket)
        self.directory = directory
        self.s3uri = 's3://arcadia-reference-datasets/OrthoFinder_output/' + self.filename

class FoldSeekOutputFile(MultiSpeciesFile):
    def __init__(self, multispeciesbiofiledocket, directory: str, filename = ''):
        super().__init__(filename = filename, multispeciesbiofiledocket = multispeciesbiofiledocket)
        self.directory = directory
        self.s3uri = 's3://arcadia-reference-datasets/FoldSeek_output/' + self.filename

class JointExcFile(MultiSpeciesFile):
    def __init__(self, multispeciesbiofiledocket, directory: str, filename = '', embedding = '', sources = []):
        super().__init__(filename = filename, multispeciesbiofiledocket = multispeciesbiofiledocket)
        self.directory = directory
        self.sources = sources
        self.embedding = embedding
        self.s3uri = 's3://arcadia-reference-datasets/' + self.embedding + '_JointExc/' + self.filename

def gxc_to_exc(sample_MSD, embedding_df, exc_file):
    import pandas as pd
    import os
    
    embedding = embedding_df.columns[0]
    
    if embedding == 'Orthogroup':
        print('Using Orthogroup embeddings as expected from OrthoFinder')
    elif embedding == 'StruCluster':
        print('Using StruCluster embeddings as expected from FoldSeek')
    else:
        raise Exception('Dataframe must have "Orthogroup" or "StruCluster" as embedding')
    
    # Iterates through all of the species in the Species BioFileDocket
    for pre in sample_MSD.species_BioFileDockets.keys():
    
        # Generates filename automatically
        exc_filename = sample_MSD.species_BioFileDockets[pre].gxc.filename.replace('.' + sample_MSD.species_BioFileDockets[pre].gxc.filetype, '_as' + embedding + '.' + sample_MSD.species_BioFileDockets[pre].gxc.filetype)
    
        # Generates file object
        exc = ExcFile(
            filename = exc_filename,
            sampledict = sample_MSD.SampleDicts[pre],
            gxcfile = sample_MSD.species_BioFileDockets[pre].gxc,
            embedding = embedding
            )
        embedding_exc = embedding + '_exc'
    
        # Checks whether an Embedding_exc file already exists; avoid re-generating if it does
        if os.path.exists(exc.path):
            print(embedding_exc + 'file already exists at', exc.path, 'skipping')
            sample_MSD.species_BioFileDockets[pre].add_keyfile(exc, embedding_exc)
        
            continue
        
        # Copies orthogroups dataframe to do transformations without modifying original
        embedding_df_copy = embedding_df.copy(deep = True)
        
        if embedding == 'Orthogroup':
            # Automatically gets the expected column name of the OrthoFinder file
            species_column = sample_MSD.species_BioFileDockets[pre].cdna.filename + '.transdecoder'
            print('Expanding column', species_column)
    
            # Expands orthogroups column for species-specific dataset
            embedding_df_copy['protein_id'] = embedding_df_copy[species_column].str.split(', ')
            embedding_df_copy = embedding_df_copy.explode('protein_id')
            embedding_df_copy['transcript_id'] = embedding_df_copy['protein_id'].str.split('.', expand = True)[0]
    
            # Gets id mapping between transcript, protein, and Orthogroup ids
            keys = embedding_df_copy[['transcript_id', 'protein_id', embedding]].drop_duplicates()
            keys.dropna(inplace = True)
            
            print('Extracting keys for', embedding)
            
            idmm_kind = 'og_idmm'
            fileouttype = 'OGfile'
            
            original_idmm = sample_MSD.species_BioFileDockets[pre].gtf_idmm
            
            # Loads the gtf id mapping matrix                                                                    
            original_idmm_df = pd.read_csv(original_idmm.path, index_col = 0, sep = '\t')
            
            # Merges original idmm with orthogroup info, generating a new idmm to be used in downstream analysis
            exc_idmm_df = original_idmm_df.merge(keys, on = 'transcript_id')  
            exc_idmm_filename = pre + '_' + sample_MSD.SampleDicts[pre].conditions + '_' + idmm_kind + '.tsv'
            
        elif embedding == 'StruCluster':
            # Automatically gets the expected column name of the FoldSeek output file
            species_column = pre
            
            # Expands struclusters column for species-specific dataset
            embedding_df_copy['uniprot_id'] = embedding_df_copy[species_column].str.split(',')
            embedding_df_copy = embedding_df_copy.explode('uniprot_id')
    
            # Gets id mapping between transcript and FoldSeek ids
            keys = embedding_df_copy[['uniprot_id', embedding]].drop_duplicates()
            keys.dropna(inplace = True)
            
            idmm_kind = 'sc_idmm'
            fileouttype = 'SCfile'
            
            original_idmm = sample_MSD.species_BioFileDockets[pre].uniprot_idmm
            
            # Loads the gtf id mapping matrix                                                                    
            original_idmm_df = pd.read_csv(original_idmm.path, index_col = 0, sep = '\t')
    
            # Merges original idmm with orthogroup info, generating a new idmm to be used in downstream analysis
            exc_idmm_df = original_idmm_df.merge(keys, on = 'uniprot_id')  
            exc_idmm_filename = pre + '_' + sample_MSD.SampleDicts[pre].conditions + '_' + idmm_kind + '.tsv'
    
        exc_idmm = IdmmFile(
            filename = exc_idmm_filename,
            sampledict = sample_MSD.SampleDicts[pre],
            kind = idmm_kind,
            sources = [exc_file, original_idmm]
            )
        print('Saving keys to', exc_idmm_filename)
        exc_idmm_df.to_csv(exc_idmm.path, sep = '\t', index = None)
        print('Done saving', exc_idmm_filename)
    
        # Adds the new og_idmm to the BioFileDocket for the species
        sample_MSD.species_BioFileDockets[pre].add_keyfile(exc_idmm, idmm_kind)
        sample_MSD.species_BioFileDockets[pre].add_keyfile(exc_file, sample_MSD.species_concat + '_' + fileouttype)
        
        # Extracts gene_name to embedding mapping keys
        gene_keys = exc_idmm_df[['gene_name', embedding]].drop_duplicates()
        
        print('Generating exc file at', exc.path)
        # Reads in original gxc matrix file
        gxc_df = pd.read_csv(sample_MSD.species_BioFileDockets[pre].gxc.path, sep = '\t')
        # Automatically gets the first column name of file for later use
        gxc_original_dataname = gxc_df.columns[0]
        # Renames that column to 'gene_name' for easier merging
        gxc_df.rename(columns = {gxc_original_dataname: 'gene_name'}, inplace = True)
        # Merges gxc with orthogroup gene keys to generate exc dataframe
        exc_df = gene_keys.merge(gxc_df, on = 'gene_name')
    
        # Removes 'gene_name' column
        exc_df = exc_df.drop(columns = 'gene_name')
        # Aggregates read counts per cell by Orthogroup ID
        exc_df = exc_df.groupby(embedding).agg({i: ('first' if i == embedding else 'sum') for i in exc_df.columns}).reset_index(drop = True)
        exc_df = exc_df.drop_duplicates(keep = 'first', subset = exc_df.columns[1:])
        
        print('Preview of exc file:', exc.path)
        display(exc_df)
    
        # Saves new exc matrix to file and puts it into the species BioFileDocket
        exc_df.to_csv(exc.path, sep = '\t', index = None)
        print('Exc file saved at', exc.path)
        
        sample_MSD.species_BioFileDockets[pre].add_keyfile(exc, embedding_exc)