from string_functions import *
import os
import subprocess

GIT_HOME = subprocess.run(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).stdout.decode("utf-8").strip('\n')
GLOBAL_OUTPUT_DIRECTORY = GIT_HOME + '/output/'

if not os.path.exists(GLOBAL_OUTPUT_DIRECTORY):
    os.mkdir(GLOBAL_OUTPUT_DIRECTORY)
    print(GLOBAL_OUTPUT_DIRECTORY, 'does not exist; making it now')

def s3_transfer(to_loc, from_loc):
    subprocess.run(['aws', 's3', 'cp', to_loc, from_loc])

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
        return 's3://arcadia-reference-datasets/glial-origins-pkl/' + self.dill_filename
    
    @property
    def dill_filename(self):
        return '_'.join([prefixify(self.species), self.conditions, 'sample_BioFileDocket.pkl'])
    
    @property
    def dill_filepath(self):
        return self.directory + self.dill_filename
    
    @property
    def sampledict(self):
        return SampleDict(self.species, self.conditions, self.directory)
       
    def set_taxid(self, taxid):
        self.taxid = str(taxid)
        
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
    
    def remove_keyfile(self, key):
        delattr(self, key)
            
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
    
    def local_to_s3(self, overwrite = False):
        files = {i:j for i,j in dict(vars(self)).items() if isinstance(j, BioFile)}

        for file in files.values():
            file.push_to_s3(overwrite)
    
    def s3_to_local(self, overwrite = False):
        files = {i:j for i,j in dict(vars(self)).items() if isinstance(j, BioFile)}
        
        for file in files.values():
            file.get_from_s3(overwrite)
    
    def get_from_s3(self, overwrite = False):
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
    def __init__(self, sampledict: SampleDict, filename = '', **kwargs):
        self.filename = filename
        self.species = sampledict.species
        self.conditions = sampledict.conditions
        self.directory = sampledict.directory
        self.s3uri = None
        
        if 'url' and 'protocol' in kwargs.keys():
            self.get_from_url(kwargs['url'], kwargs['protocol'], filename = filename)
        elif 's3uri' in kwargs.keys():
            self.get_from_s3(kwargs['s3uri'])
    
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

    def get_from_url(self, url: str, protocol: str, filename = ''):
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
        output = GenomeGtfFile(filename, sampledict = self.sampledict, GenomeFastaFile = self.reference_genome)
        
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