# utils/biofile_handling.py

"""BioFile objects and functions for managing collections of biological data."""

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

class metadata_object(object):
    """Simple dummy object to enable dot access to keys."""
    def add(self, key: str, value: any, replace = True):
        """Adds a metadata feature using a unique key identifier.

        >- Checks to make sure that key does not already exist.  
        >- Also makes sure key is only alphanumeric or underscores.  
        
        Args:
            key (str): a unique key
            value (obj): any field you want to record
        """
        if not re.match(r'^\w+$', key):
            raise Exception('key can only include alphanumeric characters and underscores')
        if hasattr(self, key):
            if replace:
                print('overwriting ' + key)
            else:
                print('key "' + key + '" already exists, ignoring')
                return
        setattr(self, key, value)

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
        metadata (metadata_object): a collector for miscellaneous metadata.
        files (dict): a dictionary of assorted files associated with the data.
        keyfiles (attr): unique key-identified attributes representing key files.
            You can add these using BioFileDocket.add_keyfile() or .add_keyfiles().
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
        self.metadata = metadata_object()
       
    @property
    def species_prefix(self):
        """`str`: runs prefixify(species)."""
        return prefixify(self.species)
    
    @property
    def s3uri(self):
        """`str`: the S3 URI of the BioFileDocket."""
        return 's3://arcadia-reference-datasets/glial-origins-pkl/' + self.dill_filename
    
    @property
    def dill_filename(self):
        """`str`: the unique filename of the BioFileDocket .pkl file."""
        return '_'.join([prefixify(self.species), self.conditions, 'BioFileDocket.pkl'])
    
    @property
    def dill_filepath(self):
        """`str`: the path of the BioFileDocket .pkl file."""
        return self.directory + self.dill_filename
    
    @property
    def sampledict(self):
        """`SampleDict`: the full SampleDict of the BioFileDocket."""
        return SampleDict(self.species, self.conditions, self.directory)
       
    def set_taxid(self, taxid: Union[str, int]):
        """Adds a taxid attribute to the BioFileDocket.  
        
        Tolerates either `int` or `str` input.
        """
        self.metadata.add('taxid', str(taxid))
        
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
            
    def add_keyfile(self, key: str, BioFile, overwrite = False):
        """Adds a BioFile object using a unique key identifier.  
        
        Checks to make sure that key does not already exist.
        Also makes sure key is only alphanumeric or underscores.
        
        Args:
            key (str): a unique key
            BioFile (BioFile): a BioFile object
        """
        if not re.match(r'^\w+$', key):
            raise Exception('key can only include alphanumeric characters and underscores')
        if hasattr(self, key) and not overwrite:
            print('key "' + key + '" already exists, ignoring')
        else:
            setattr(self, key, BioFile)
        
    def add_keyfiles(self, dictionary: dict, overwrite = False):
        """Adds a dictionary of BioFile objects using key:BioFile object pairs."""
        for key in dictionary:
            self.add_keyfile(key, dictionary[key], overwrite)
    
    def remove_keyfile(self, key, warn = True):
        """Deletes a keyfile from the BioFileDocket if warn == `False`."""
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
            self (BioFileDocket): returns the unpickled BioFileDocket object.
        
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

class MultiSpeciesBioFileDocket(BioFileDocket):
    """MultiSpeciesBioFileDocket objects collect BioFileDockets and relate them to one another.
    
    Args:
        species_dict (dict): key is species name in 'Genus_species' format; value is conditions.
            These must be exact matches, or it will not work.
        global_conditions (str): a summary identifier for this collection of species' datasets.
        analysis_type (str): one-word description of the analysis type.
    
    Attributes:
        directory (str): output directory for files in the dataset.
    """
    def __init__(self, species_dict: dict, global_conditions: str, analysis_type: str):
        self.species_dict = species_dict
        self.global_conditions = global_conditions
        self.analysis_type = analysis_type
        
        # Dummy variables to be backwards-compatible with BioFileDocket
        self.species = self.species_concat
        self.conditions = '_'.join([self.global_conditions, self.analysis_type])
        
        self.metadata = metadata_object()
        
        directory = make_output_directory(self.species_concat, self.conditions)
        self.directory = directory
    
    @property
    def species_concat(self):
        """`str`: concatenation of species prefixes in alphabetical order."""
        return ''.join(sorted([prefixify(species) for species in self.species_dict]))
    
    @property
    def sampledict(self):
        """`SampleDict`: a MultiSpeciesBioFileDocket-formatted sampledict."""
        return SampleDict(self.species_concat, self.conditions, self.directory)
    
    @property
    def dill_filename(self):
        """`str`: the unique filename of the BioFileDocket .pkl file."""
        return '_'.join([self.species_concat, self.global_conditions, self.analysis_type, 'MultiSpeciesBioFileDocket.pkl'])
    
    def get_BioFileDockets(self):
        """Gets species BioFileDocket .pkl files from S3 for each species in the species_dict."""
        import dill
        
        self.species_BioFileDockets = {}
        
        for species in self.species_dicts:
            conditions = self.species_dicts[species]
            species_prefix = prefixify(species)
            species_directory = make_output_directory(species, conditions, stringonly = True)
            
            dill_filename = '_'.join([species_prefix, conditions, 'BioFileDocket.pkl'])
            dill_filepath = species_directory + dill_filename
            dill_s3uri = 's3://arcadia-reference-datasets/glial-origins-pkl/' + dill_filename
            
            if not os.path.exists(dill_filepath):
                s3_transfer(dill_s3uri, dill_filepath)
            
            with open(dill_filepath, 'rb') as file:
                self.species_BioFileDockets[species_prefix] = dill.load(file)
                
        return self
    
    def s3_to_local(self, overwrite = False):
        """Iteratively calls s3_to_local on all BioFileDockets in the group."""
        for species_prefix in self.species_BioFileDockets:
            self.species_BioFileDockets[species_prefix].s3_to_local(overwrite)
    
    def local_to_s3(self, overwrite = False):
        """Iteratively calls local_to_s3 on all BioFileDockets in the group."""
        for species_prefix in self.species_BioFileDockets:
            self.species_BioFileDockets[species_prefix].local_to_s3(overwrite)
        
# Class BioFile is used to carry metadata for each file
class BioFile:
    """BioFile objects collect metadata about biological filetypes.
    
    Args:
        sampledict (SampleDict): a SampleDict object from the BioFileDocket.
        filename (str, optional): the name of the file.
        url (str, optional): when downloading a file on object creation, pass a string url along with a protocol.
        protocol (str, optional): passed along with a url for automatic download on object creation.
        s3uri (str, optional): the s3uri of the file, if downloading from s3 upon object creation.
        unzip (bool, optional): whether or not to unzip a file on download. Defaults to True.
    """
    def __init__(self, sampledict: SampleDict, filename = '', url = None, protocol = None, s3uri = None, unzip = True):
        self.filename = filename
        self.species = sampledict.species
        self.conditions = sampledict.conditions
        self.directory = sampledict.directory
        self.s3uri = s3uri
        self.metadata = metadata_object()
        
        if url != None and protocol != None:
            self.get_from_url(url = url, protocol = protocol, filename = filename, unzip = unzip)
        elif s3uri != None:
            self.get_from_s3(s3uri)
        elif filename == '':
            raise Exception('File must minimally have a filename.\nTo download from url, pass both a url and protocol.\nTo download from s3uri, pass s3uri.')
            
    @property
    def exists(self):
        """`bool`: checks whether the file currently exists."""
        return os.path.exists(self.path)
    
    @property
    def path(self):
        """`str`: path to the file, including the filename."""
        return self.directory + self.filename
    
    @property
    def species_prefix(self):
        """`str`: runs `prefixify(species)`."""
        return prefixify(self.species)
    
    @property
    def sampledict(self):
        """`SampleDict`: a SampleDict object for the file."""
        return SampleDict(self.species, self.conditions, self.directory)
    
    @property
    def filetype(self):
        """`str`: infers filetype using the string after the final `.` in filename."""
        return self.filename.split('.')[-1]
    
    def add_s3uri(self, s3uri: str):
        """Adds an s3uri to the BioFile object if it doesn't already exist."""
        if self.s3uri == None:
            self.s3uri = s3uri
        else:
            raise Exception('This file already has an S3 URI at ' + self.s3uri)
        return None
    
    def get_from_s3(self, overwrite = False):
        """Downloads the BioFile from AWS S3.
        
        Args:
            overwrite (bool): decide whether to overwrite existing files. Defaults to False.
        """
        import os
        import subprocess
        
        if self.filename == '':
            self.filename = self.s3uri.split('/')[-1]
            print('inferring file name as', self.filename)
        
        if not self.exists:
            s3_transfer(self.s3uri, self.path)
            return self
        elif overwrite:
            s3_transfer(self.s3uri, self.path)
            return self
        else:
            print('file', self.filename, 'already exists at', self.path)
            return self
    
    def push_to_s3(self, overwrite = False):
        """Uploads the BioFile to AWS S3.
        
        Args:
            overwrite (bool): decide whether to overwrite existing files. Defaults to False.
        """
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
        """Downloads the BioFile from a URL using a chosen protocol, unzipping optionally.
        
        Args:
            url (str): url of the file.
            protocol (str): protocol to be used for download.
            filename (str, optional): name of file to be saved. If empty, will generate name from URL.
            unzip (bool): decide whether to unzip if it is a zipped file. Defaults to True.
        """
        import os, subprocess
        
        self.url = url

        protocols = ['curl', 'wget']
        if protocol not in protocols:
            raise ValueError("Invalid protocol. Expected one of: %s" % protocols)   
            
        if filename == '':
            filename = self.url.split('/')[-1]
            print('inferring file name as', filename)
            self.filename = filename
            
        if self.exists:
            print('file', self.filename, 'already exists at', self.path)
            
            if unzip:
                self.unzip()
                
            return self
        else:
            if protocol == 'curl':
                subprocess.run([protocol, url, '--output', self.path])
            elif protocol == 'wget':
                subprocess.run([protocol, '-O', self.path, self.url])
            print('downloaded file', self.path)
            
            if unzip:
                self.unzip()
            
            return self
        
    def unzip(self):
        """Unzips files ending in .gz, .gzip, or .zip."""
        compressions = ['gz', 'gzip', 'zip']
        compression = self.filename.split('.')[-1]
        
        if compression not in compressions:
            return self
        else:
            output_name = self.filename.replace('.' + compression, '')
            output_loc = self.directory + output_name
            
        if compression == 'gz' or 'gzip':
            subprocess.run(['gunzip', output_loc])
        if compression == 'zip':
            subprocess.run(['unzip', output_loc])
        
        print('file', self.filename, 'unzipped and object renamed to', output_name)
        self.filename = output_name
        
        return self

class GenomeFastaFile(BioFile):
    """A BioFile object for genome fasta files.
    
    Args:
        version (str): the specific version number of the genome.
    """
    def __init__(self, version: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.version = version
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
    
    def rename_RefSeq_chromosomes(self, replace = False):
        """Renames the chromosomes of a RefSeq genome to numerical chromosomes."""
        from Bio import SeqIO

        new_filename = self.filename + '.renamed.fa'
        new_filepath = self.directory + new_filename
        
        new_file = GenomeFastaFile(
                filename = new_filename,
                sampledict = self.sampledict,
                version = self.version
            )
        if new_file.exists:
            return new_file

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
            return new_file
        
    def get_transdecoder_cdna_gtf(self, GenomeGtfFile, TRANSDECODER_LOC, **kwargs):
        """Generates a TransdecoderCdnaFile given a GenomeGtfFile.
        
        Args:
            GenomeGtfFile (GenomeGtfFile): a GenomeGtfFile object that is compatible with the GenomeFastaFile.
            TRANSDECODER_LOC (str): the path to the parent directory of Transdecoder.
        
        Returns:
            cdna_output (TransdecoderCdnaFile): a TransdecoderCdnaFile object.
        """
        PERL_SCRIPT_LOC = TRANSDECODER_LOC + 'util/gtf_genome_to_cdna_fasta.pl'
        
        import os
        import subprocess
        
        cdna_name = self.filename.replace('.' + self.filetype, '_cDNA.' + self.filetype)
        cdna_output = TransdecoderCdnaFile(
            sampledict = self.sampledict,
            filename = cdna_name,
            reference_genome = self,
            reference_annot = GenomeGtfFile
        )
        if cdna_output.exists:
            return cdna_output
        else:
            output_file = open(cdna_output.path, "w")
            
            subprocess.call([PERL_SCRIPT_LOC, GenomeGtfFile.path, self.path], stdout=output_file)
            return cdna_output

class TransdecoderCdnaFile(BioFile):
    """A BioFile object for Transdecoder cDNA files.
    
    Args:
        reference_genome (GenomeFastaFile): the BioFile object of the associated genome fasta file.
        reference_annot (GenomeGffFile or GenomeGtfFile): the BioFile object of the associated genome annotation file.
    """
    def __init__(self, reference_genome, reference_annot, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/transcriptome/' + self.filename
        self.reference_genome = reference_genome
        self.reference_annot = reference_annot

    def to_pep_files(self, TDLONGORF_LOC, TDPREDICT_LOC):
        """ Generates a peptide file from a Transdecoder cDNA file.
        
        Args:
            TDLONGORF_LOC (str): path to the TransDecoder.LongOrfs binary.
            TDPREDICT_LOC (str): path to the TransDecoder.Predict binary.
        
        Returns:
            output_dict (dict of TransdecoderOutFile): a dictionary of Transdecoder output files.
        """
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
                reference_genome = self.reference_genome, 
                reference_annot = self.reference_annot,
                reference_cDNA = self
            )
            subprocess.run(['mv', output_file.filename, output_file.path])
            output_dict['transdecoder_' + suffix.replace('.', '')] = output_file
        
        return output_dict
    
class TransdecoderOutFile(BioFile):
    """A BioFile object for Transdecoder cDNA files.
    
    Args:
        reference_genome (GenomeFastaFile): the BioFile object of the associated genome fasta file.
        reference_annot (GenomeGffFile or GenomeGtfFile): the BioFile object of the associated genome annotation file.
        reference_cDNA (TransdecoderCdnaFile): the BioFile object of the associated Transdecoder cDNA file.
    """
    def __init__(self, reference_genome, reference_annot, reference_cDNA, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/proteome/' + self.filename
        self.reference_genome = reference_genome
        self.reference_annot = reference_annot
        self.reference_cDNA = reference_cDNA

class GenomeGffFile(BioFile):
    """A BioFile object for a genome reference in GFF format.
    
    If a file ends in `.gff3`, renames it to end in `.gff` instead.
    
    Args:
        reference_genome (GenomeFastaFile): the BioFile object of the associated genome fasta file.
    """
    def __init__(self, reference_genome, *args, **kwargs):
        import subprocess, os
        
        super().__init__(*args, **kwargs)
        
        if self.filename.split('.')[-1] == 'gff3':
            new_filename = self.filename.replace('gff3', 'gff')
            new_path = self.directory + new_filename
            subprocess.run(['mv', self.path, new_path])
            
            print('Renaming GFF3 file', self.filename, 'to', new_filename)
            self.filename = new_filename
        
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = reference_genome
    
    def to_gtf(self, GFFREAD_LOC, keep_all = False):
        import subprocess
        import os
        
        filename = self.filename.replace('.gff', '.gtf')
        output = GenomeGtfFile(filename = filename, sampledict = self.sampledict, reference_genome = self.reference_genome)
        
        if output.exists:
            print('Converted file', output.filename, 'already exists at:\n', output.path)
            return output
        
        subprocess.run([GFFREAD_LOC, self.path, '-T', '-o', output.path])
        return output

class GenomeGtfFile(BioFile):
    """A BioFile object for a genome reference in GTF format.
    
    Args:
        reference_genome (GenomeFastaFile): the BioFile object of the associated genome fasta file.
    """
    def __init__(self, reference_genome, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = reference_genome

class CellRangerFileGroup:
    """A special object that collects the multiple files of a CellRanger file.
    
    Args:
        sampledict (SampleDict): a SampleDict object from the BioFileDocket.
        barcodes_address (str): the address of the barcodes file as either url or s3uri.
        features_address (str): the address of the features file as either url or s3uri.
        matrix_address (str): the address of the matrix file as either url or s3uri.
        how (str): 'url' or 's3uri' to choose a method of download.
        protocol (str): when using url, the protocol to use. Defaults to 'curl'.
    """
    def __init__(self, sampledict: SampleDict, barcodes_address: str, features_address: str, matrix_address: str, how = None, protocol = 'curl'):
        self.sampledict = sampledict
        
        if how == 'url':
            self.barcodes = CellRangerBarcodesFile(sampledict = sampledict, url = barcodes_address, protocol = protocol)
            self.matrix = CellRangerMatrixFile(sampledict = sampledict, url = matrix_address, protocol = protocol)
            self.features = CellRangerFeaturesFile(sampledict = sampledict, url = features_address, protocol = protocol)
        elif how == 's3uri':
            self.barcodes = CellRangerBarcodesFile(sampledict = sampledict, s3uri = barcodes_address)
            self.matrix = CellRangerMatrixFile(sampledict = sampledict, s3uri = matrix_address)
            self.features = CellRangerFeaturesFile(sampledict = sampledict, s3uri = features_address)
        else:
            raise Exception('Must specify how = "url" or "s3uri".')

    # Converts a group of CellRanger files into a single dense csv matrix
    def to_gxc(self, reference_genome, reference_annot, filename: str, overwrite = False):
        """Converts the group of cell ranger files to a GxcFile and returns the file object."""
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
        output = GxcFile(reference_genome = reference_genome, reference_annot = reference_annot, sampledict = self.sampledict, filename = filename)
        
        if not output.exists:
            matrix.to_csv(output.path, index=False, sep = '\t')
        else:
            print('File already exists at', output.path, '. Set overwrite = True to overwrite.')
        
        return output

class CellRangerBarcodesFile(BioFile):
    """A BioFile object for the barcodes file of a CellRanger file."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        
class CellRangerMatrixFile(BioFile):
    """A BioFile object for the matrix file of a CellRanger file."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename

class CellRangerFeaturesFile(BioFile):
    """A BioFile object for the features file of a CellRanger file."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename

class LoomFile(BioFile):
    """A BioFile object for a Loom file.
    
    Not built out yet!"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
    
    def to_gxc(self, filename, reference_genome, reference_annot, overwrite = False):
        print('This is not built yet!')
        
class GxcFile(BioFile):
    """A BioFile object for a genes x cells matrix.
    
    Args:
        reference_genome (GenomeFastaFile): the BioFile object of the associated genome fasta file.
        reference_annot (GenomeGffFile or GenomeGtfFile): the BioFile object of the associated genome annotation file.
    """
    def __init__(self, reference_genome, reference_annot, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.reference_genome = reference_genome
        self.reference_annot = reference_annot

class ExcFile(BioFile):
    """A BioFile object for an embedding x cells matrix.
    
    Args:
        gxcfile (GxcFile): the original `GxcFile` object used to generate the new embedding.
        embedding (str): a description of the new embedding.
    """
    def __init__(self, gxcfile, embedding: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.original = gxcfile
        self.embedding = embedding

class IdmmFile(BioFile):
    """An ID-mapping matrix file. Each row represents a feature and each column represents a name of that feature in a namespace.
    
    Args:
        sources (list of BioFile objects): a list of the source BioFile objects.
        kind (str): a one-word description of the kind of idmm.
    """
    def __init__(self, kind: str, sources: list, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        self.kind = kind
        
class CellAnnotFile(BioFile):
    """A cell annotation matrix file with two columns: 'cell_barcode' and 'celltype'.
    
    Each cell should have a cell barcode exactly matching its name in the Gxc or Exc embedding.  
    
    Args:
        sources (list of BioFile objects): a list of the source BioFile objects.
    """
    def __init__(self, sources: list, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.sources = sources

class UniprotIDMapperFile(IdmmFile):
    """An BioFile object for the output file generated by calling the UNIPROT ID mapping API.
    
    Args:
        from_type (str): the source datatype (e.g. 'Ensembl')
        to_type (str): the destination datatype (e.g. 'UniprotKB')
    """
    def __init__(self, from_type: str, to_type: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.from_type = from_type
        self.to_type = to_type

class UniProtTaxidListFile(BioFile):
    """An BioFile object for pulling files from Uniprot based on input taxid.
    
    Defaults to creating a file based on a taxid passed on object creation.  
    
    Args:
        taxid (str or int): the taxid of the species of interest.
        make (bool): whether or not to make the file upon object creation. defaults to True.
    """
    def __init__(self, taxid: Union[str, int], make = True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.taxid = str(taxid)
        
        if filename == '':
            self.filename = '_'.join([self.species_prefix, 'uniprot-taxid', self.taxid, 'genes.tsv'])
        
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        
        if make:
            self.get_proteins()
        else:
            print('make set to False; skipping file creation. Use .get_proteins() to create the file separately.')
        
    def get_proteins(self):
        """Queries the Uniprot API to extract all proteins for the taxid associated with this object."""
        if not self.exists:
            with open(self.path, "w") as outfile:
                subprocess.run(['curl', 
                  'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_id%2Corganism_name&format=tsv&query=%28%28taxonomy_id%3A' + self.taxid + '%29%29'],
                  stdout = outfile)
        return self
        
class GeneListFile(BioFile):
    """An BioFile object for a gene list file, which is a .txt file with a gene on each line.
    
    Defaults to creating a file based on a list of genes passed on object creation.  
    
    Args:
        sources (list of BioFile objects): a list of the source BioFile objects.
        identifier (str): a description of the identifier type for this gene list.
        make (bool): whether or not to make the file upon object creation. defaults to True.
    """
    def __init__(self, sources: list, genes: list, identifier: str, filename = '', autoname = True, make = True, *args, **kwargs):
        if filename == '':
            filename = 'placeholder.txt'
        super().__init__(filename = filename, *args, **kwargs)
        
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        self.identifier = identifier
        
        if autoname:
            self.filename = '_'.join([prefixify(self.sampledict.species), self.sampledict.conditions, self.identifier, 'ids.txt'])
        
        if make:
            self.make_file(genes)
        else:
            print('make set to False; skipping file creation. Use .make_file(genes) to create the file separately.')
    
    def make_file(self, genes: list):
        """Makes a gene list file at a given location when passed a list object."""
        make_gene_list(genes, self.path)
    
    def get_uniprot_ids(self, ID_MAPPER_LOC: str, from_type: str, to_type: str):
        """Queries the Uniprot API to get ID mapping based on from_type and to_type, returning an UniprotIDMapperFile object.
        
        If the number of genes in the file exceeds 100k, generates temporary files to query the API multiple times.
        
        Args:
            ID_MAPPER_LOC (str): path to the ID mapping script, usually at ID_MAPPER_LOC.
            from_type (str): the source datatype (e.g. 'Ensembl')
            to_type (str): the destination datatype (e.g. 'UniprotKB')
        """
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
        """Don't use this, it's not built out yet."""
        import subprocess
        
        filename = self.filename.replace('_ids.txt', '_UniProtIDs.txt')
        
        subprocess.run([ID_MAPPER_LOC, self.path, from_type, to_type, taxid])
        output = UniprotIDMapperFile(
            filename = filename, sampledict = self.sampledict, kind = 'UniprotIDMapper', 
            sources = self.sources, from_type = from_type, to_type = to_type)
        
        return output

class MultiSpeciesFile(BioFile):
    """A special class of BioFile objects for files with multiple associated species.
    
    This subclass uses `species_concat` as its species and its conditions is a '_' concatenation of its `global_conditions` and `analysis_type`.
    
    Args:
        species_dict (dict): species dictionary from the MultiSpeciesBioFileDocket.
    """
    def __init__(self, species_dict: dict, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.species_dict = species_dict
        self.global_conditions = sampledict.conditions.split('_')[0]
        self.analysis_type = sampledict.conditions.split('_')[1]
    
    @property
    def species_concat(self):
        """str: concatenation of species prefixes in alphabetical order."""
        return ''.join(sorted([prefixify(species) for species in self.species_dict]))
    
    @property
    def sampledict(self):
        """SampleDict: a MultiSpeciesBioFileDocket-formatted sampledict."""
        return SampleDict(self.species_concat, self.conditions, self.directory)
        
class OrthoFinderOutputFile(MultiSpeciesFile):
    """A MultiSpeciesFile for OrthoFinder output."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/OrthoFinder_output/' + self.filename

class FoldSeekOutputFile(MultiSpeciesFile):
    """A MultiSpeciesFile for FoldSeek output."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/FoldSeek_output/' + self.filename

class JointExcFile(MultiSpeciesFile):
    """A MultiSpeciesFile for collecting an ExcFile with multiple species.
    
    Args:
        sources (list of GxcFile): the original `GxcFile` object used to generate the new embedding.
        embedding (str): a description of the new embedding.
    """
    def __init__(self, sources: list, embedding: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sources = sources
        self.embedding = embedding
        self.s3uri = 's3://arcadia-reference-datasets/' + self.embedding + '_JointExc/' + self.filename

def gxc_to_exc(sample_MSD, embedding_df, exc_file):
    """Converts an GxcFile into an ExcFile, returning the new file object."""
    
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
        if exc.exists:
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