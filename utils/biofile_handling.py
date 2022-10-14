from string_functions import *
import os
import subprocess

GIT_HOME = subprocess.run(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).stdout.decode("utf-8").strip('\n')
GLOBAL_OUTPUT_DIRECTORY = GIT_HOME + '/output/'
S3_BUCKET = 's3://arcadia-reference-datasets'
S3_PKL_DIR = 'glial-origins-pkl'

if not os.path.exists(GLOBAL_OUTPUT_DIRECTORY):
    os.mkdir(GLOBAL_OUTPUT_DIRECTORY)
    print(GLOBAL_OUTPUT_DIRECTORY, 'does not exist; making it now')

# Class definitions to handle BioFiles between scripts
class SampleDict:
    def __init__(self, species: str, conditions: str, directory: str):
        self.species = species
        self.conditions = conditions
        self.directory = directory

class BioFileDocket:                             
    def __init__(self, sampledict):
        self.species = sampledict.species
        self.conditions = sampledict.conditions
        self.directory = sampledict.directory
        self.files = dict()
    
    @property
    def s3uri(self):
        return '/'.join([S3_BUCKET, S3_PKL_DIR, self.dill_filename])
    
    @property
    def dill_filename(self):
        return '_'.join([prefixify(self.species), self.conditions, 'sample_BioFileDocket.pkl'])
    
    @property
    def dill_filepath(self):
        return self.directory + self.dill_filename
        
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
    
    def pickle(self):
        import dill

        with open(self.dill_filepath, 'wb') as file:
            dill.dump(self, file)
    
    def unpickle(self):
        import dill
        import os
        
        if not os.path.exists(self.dill_filepath):
            raise Exception("Can't unpickle file; .pkl file doesn't exist yet at " + self.dill_filepath)

        with open(self.dill_filepath, 'rb') as file:
            self = dill.load(file)
        
        return self
    
    def s3_sync(self, how: str, overwrite = False):
        files = {i:j for i,j in dict(vars(self)).items() if isinstance(j, BioFile)}
        
        if how == 'push':
            for file in files.values():
                file.push_to_s3(overwrite)
        elif how == 'pull':
            for file in files.values():
                file.get_from_s3(overwrite)
        else:
            raise Exception('"how=" expects "push" or "pull" as input value')
    
    def get_from_s3(self, overwrite = False):
        import os
        import subprocess
        
        if not os.path.exists(self.dill_filepath):
            subprocess.run(['aws', 's3', 'cp', self.s3uri, self.dill_filepath])
            return self
        elif overwrite:
            subprocess.run(['aws', 's3', 'cp', self.s3uri, self.dill_filepath])
            return self
        else:
            print('file', self.dill_filename, 'already exists at', self.dill_filepath)
            return self
    
    def push_to_s3(self, overwrite = False):
        import subprocess
        
        bucket = self.s3uri.lstrip('s3://').split('/')[0]
        file_key = self.s3uri.split(bucket)[1].lstrip('/')
        
        # check if file exists in S3 already
        output = subprocess.run(['aws', 's3api', 'head-object', '--bucket', bucket, '--key', file_key], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        if 'Not Found' in str(output.stderr):
            subprocess.run(['aws', 's3', 'cp', self.dill_filepath, self.s3uri])
        elif overwrite:
            print(self.dill_filename, 'already exists in S3 bucket; overwriting.')
            subprocess.run(['aws', 's3', 'cp', self.dill_filepath, self.s3uri])
        else:
            print(self.dill_filename, 'already exists in S3 bucket, skipping upload. set overwrite = True to overwrite the existing file.')
        return None

class MultiSpeciesBioFileDocket:
    def __init__(self, species_dict: dict, global_conditions: str, analysis_type):
        self.species_dict = species_dict
        self.global_conditions = global_conditions
        self.analysis_type = analysis_type
        self.SampleDicts = {}
        
        for species in species_dict:
            species_prefix = prefixify(species)
            conditions = species_dict[species]
            output_folder = '../../output/' + species_prefix + '_' + conditions + '/'
            species_SampleDict = SampleDict(species, conditions, output_folder)
            self.SampleDicts[species_prefix] = species_SampleDict
        
        self.species_concat = ''.join(sorted(self.SampleDicts.keys()))
        self.directory = '../../output/' + self.species_concat + '_' + self.global_conditions + '_' + self.analysis_type + '/'
    
    def make_directory(self):
        import os
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        return None
    
    def get_BioFileDockets(self):
        import dill
        
        self.species_BioFileDockets = {}
        
        for sampledict in self.SampleDicts.values():
            species_prefix = prefixify(sampledict.species)
            
            dill_filename = '_'.join([species_prefix, sampledict.conditions, 'sample_BioFileDocket.pkl'])
            dill_filepath = sampledict.directory + dill_filename
            dill_s3uri = 's3://arcadia-reference-datasets/glial-origins-pkl/' + dill_filename
            
            if not os.path.exists(dill_filepath):
                subprocess.run(['aws', 's3', 'cp', dill_s3uri, dill_filepath])
            
            with open(dill_filepath, 'rb') as file:
                self.species_BioFileDockets[species_prefix] = dill.load(file)
        return None
    
    def s3_sync(self, how: str, overwrite = False):
        for species_prefix in self.species_BioFileDockets:
            self.species_BioFileDockets[species_prefix].s3_sync(how = how, overwrite = overwrite)
        
        
# Class BioFile is used to carry metadata for each file
class BioFile:
    def __init__(self, filename: str, sampledict: SampleDict, **kwargs):
        self.filename = filename
        self.species = sampledict.species
        self.conditions = sampledict.conditions
        self.directory = sampledict.directory
        self.sampledict = SampleDict(self.species, self.conditions, self.directory)
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
    
    def get_from_s3(self, overwrite = False):
        import os
        import subprocess
        
        if not os.path.exists(self.path):
            subprocess.run(['aws', 's3', 'cp', self.s3uri, self.path])
            return self
        elif overwrite:
            subprocess.run(['aws', 's3', 'cp', self.s3uri, self.path])
            return self
        else:
            print('file', self.filename, 'already exists at', self.path)
            return self
    
    def push_to_s3(self, overwrite = False):
        import subprocess
        
        if self.s3uri is None:
            raise Exception('There is no s3uri. Use add_s3uri() to add a URI.')
        
        bucket = self.s3uri.lstrip('s3://').split('/')[0]
        file_key = self.s3uri.split(bucket)[1].lstrip('/')
        
        # check if file exists in S3 already
        # suppress stdout and save stderr as file to detect non-existing files
        output = subprocess.run(['aws', 's3api', 'head-object', '--bucket', bucket, '--key', file_key], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        if 'Not Found' in str(output.stderr):
            subprocess.run(['aws', 's3', 'cp', self.path, self.s3uri])
        elif overwrite:
            print(self.filename, 'already exists in S3 bucket; overwriting.')
            subprocess.run(['aws', 's3', 'cp', self.path, self.s3uri])
        else:
            print(self.filename, 'already exists in S3 bucket, skipping upload. set overwrite = True to overwrite the existing file.')
        
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
            
            if os.path.exists(self.path):
                return self
            
            if protocol == 'curl':
                subprocess.run([protocol, url, '--output', self.path])
            return self

class GenomeFastaFile(BioFile):
    def __init__(self, filename: str, sampledict: SampleDict, version: str, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.version = version
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/genome/' + self.filename
        
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
    def __init__(self, filename: str, sampledict: SampleDict, GenomeFastaFile, GenomeAnnotFile, **kwargs):
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
                output_filename, 
                self.sampledict, 
                GenomeFastaFile = self.reference_genome, 
                GenomeAnnotFile = self.reference_annot,
                TransdecoderCdnaFile = self
            )
            subprocess.run(['mv', output_file.filename, output_file.path])
            output_dict['transdecoder_' + suffix.replace('.', '')] = output_file
        
        return output_dict
    
class TransdecoderOutFile(BioFile):
    def __init__(self, filename, sampledict, GenomeFastaFile, GenomeAnnotFile, TransdecoderCdnaFile, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/proteome/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile
        self.reference_cDNA = TransdecoderCdnaFile

class GenomeGffFile(BioFile):
    def __init__(self, filename: str, sampledict: SampleDict, GenomeFastaFile: GenomeFastaFile, **kwargs):
        import subprocess
        import os
        
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        
        if self.filename.split('.')[-1] == 'gff3':
            new_filename = self.filename.replace('gff3', 'gff')
            new_path = self.directory + new_filename
            subprocess.run(['mv', self.path, new_path])
            self.filename = new_filename
            self.path = new_path
        
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = GenomeFastaFile
    
    def to_gtf(self, GFFREAD_LOC, keep_all = False):
        import subprocess
        import os
        
        filename = self.filename.replace('.gff', '.gtf')
        output = GenomeGtfFile(filename, sampledict = self.sampledict, GenomeFastaFile = self.reference_genome)
        
        if os.path.exists(output.path):
            print('Converted file', output.filename, 'already exists at:\n', output.path)
            return output
        
        subprocess.run([GFFREAD_LOC, self.path, '-T', '-o', output.path])
        
        return output

class GenomeGtfFile(BioFile):
    def __init__(self, filename, sampledict, GenomeFastaFile, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/annotation/' + self.filename
        self.reference_genome = GenomeFastaFile

class GxcFile(BioFile):
    def __init__(self, filename: str, sampledict: SampleDict, GenomeFastaFile, GenomeAnnotFile, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.reference_genome = GenomeFastaFile
        self.reference_annot = GenomeAnnotFile

class ExcFile(BioFile):
    def __init__(self, filename: str, sampledict: SampleDict, gxcfile, embedding, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.original = gxcfile
        self.embedding = embedding

class IdmmFile(BioFile):
    def __init__(self, filename, sampledict, kind, sources: list, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/genomics_reference/mapping_file/' + self.filename
        self.sources = sources
        self.kind = kind
        
class CellAnnotFile(BioFile):
    def __init__(self, filename, sampledict, sources: list, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, **kwargs)
        self.s3uri = 's3://arcadia-reference-datasets/organisms/' + self.species + '/functional_sequencing/scRNA-Seq/' + self.filename
        self.sources = sources

class UniprotIDMapperFile(IdmmFile):
    def __init__(self, filename, sampledict, kind, sources: list, from_type, to_type, **kwargs):
        super().__init__(filename = filename, sampledict = sampledict, kind = kind, sources = sources, **kwargs)
        self.from_type = from_type
        self.to_type = to_type
        
class GeneListFile(BioFile):
    def __init__(self, filename, sampledict, sources: list, genes: list, identifier: str, **kwargs):
        from string_functions import make_gene_list, prefixify
        
        self.identifier = identifier
        new_filename = '_'.join([prefixify(sampledict.species), sampledict.conditions, self.identifier, 'ids.txt'])
        
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
        
        filename = self.filename.replace('_ids.txt', '_UniProtIDs.txt')
        
        subprocess.run([ID_MAPPER_LOC, self.path, from_type, to_type])
        output = UniprotIDMapperFile(
            filename, sampledict = self.sampledict, kind = 'UniprotIDMapper', 
            sources = self.sources, from_type = from_type, to_type = to_type)
        
        return output
    
    def get_uniprot_ids_by_genename(self, ID_MAPPER_LOC, from_type, to_type, taxid):
        import subprocess
        
        filename = self.filename.replace('_ids.txt', '_UniProtIDs.txt')
        
        subprocess.run([ID_MAPPER_LOC, self.path, from_type, to_type, taxid])
        output = UniprotIDMapperFile(
            filename, sampledict = self.sampledict, kind = 'UniprotIDMapper', 
            sources = self.sources, from_type = from_type, to_type = to_type)
        
        return output

class MultiSpeciesFile(BioFile):
    def __init__(self, filename, multispeciesbiofiledocket):
        self.filename = filename
        self.multispeciesbiofiledocket = multispeciesbiofiledocket
        self.species = list(self.multispeciesbiofiledocket.species_dict.keys())
        self.global_conditions = self.multispeciesbiofiledocket.global_conditions
        self.directory = self.multispeciesbiofiledocket.directory
        self.path = self.directory + self.filename
        self.species_concat = self.multispeciesbiofiledocket.species_concat
        self.s3uri = None

class OrthoFinderOutputFile(MultiSpeciesFile):
    def __init__(self, filename, multispeciesbiofiledocket, directory: str):
        super().__init__(filename = filename, multispeciesbiofiledocket = multispeciesbiofiledocket)
        self.directory = directory
        self.path = self.directory + self.filename
        self.s3uri = '/'.join([S3_BUCKET, 'orthofinder', self.filename])