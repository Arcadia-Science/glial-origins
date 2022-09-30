from string_functions import *
import os
import subprocess

GIT_HOME = subprocess.run(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).stdout.decode("utf-8").strip('\n')
GLOBAL_OUTPUT_DIRECTORY = GIT_HOME + '/output/'

if not os.path.exists(GLOBAL_OUTPUT_DIRECTORY):
    os.mkdir(GLOBAL_OUTPUT_DIRECTORY)
    print(GLOBAL_OUTPUT_DIRECTORY, 'does not exist; making it now')

# Class definitions to handle BioFiles between scripts
class SampleDict:
    def __init__(self, species: str, conditions: str, directory: str):
        self.species = species
        self.conditions = conditions
        self.directory = directory

class Docket:                             
    def __init__(self, sampledict):
        self.species = sampledict.species
        self.conditions = sampledict.conditions
        self.directory = sampledict.directory
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
    
    def pickle(self):
        import dill
        self.dill_filepath = self.directory + '_'.join([prefixify(self.species), self.conditions, 'sample_Docket.pkl'])

        with open(self.dill_filepath, 'wb') as file:
            dill.dump(self, file)
    
    def unpickle(self):
        import dill
        import os
        
        dill_filepath = self.directory + '_'.join([prefixify(self.species), self.conditions, 'sample_Docket.pkl'])
        if os.path.exists(dill_filepath):
            self.dill_filepath = dill_filepath

        with open(self.dill_filepath, 'rb') as file:
            self = dill.load(file)
        
        return self

class MultiSpeciesDocket:
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
    
    def get_Dockets(self):
        import dill
        
        self.species_Dockets = {}
        
        for sampledict in self.SampleDicts.values():
            species_prefix = prefixify(sampledict.species)
            with open(sampledict.directory + '_'.join([species_prefix, sampledict.conditions, 'sample_Docket.pkl']), 'rb') as file:
                self.species_Dockets[species_prefix] = dill.load(file)
        return None
        
        
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
        
        if filename != '':
            print('filename is ignored and generated by input as ', new_filename)

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

class MultiSpeciesFile():
    def __init__(self, filename, multispeciesdocket):
        self.filename = filename
        self.multispeciesdocket = multispeciesdocket
        self.species = list(self.multispeciesdocket.species_dict.keys())
        self.global_conditions = self.multispeciesdocket.global_conditions
        self.directory = self.multispeciesdocket.directory
        self.path = self.directory + self.filename
        self.species_concat = self.multispeciesdocket.species_concat
        self.s3uri = None

class OrthoFinderOutputFile(MultiSpeciesFile):
    def __init__(self, filename, multispeciesdocket, directory: str):
        super().__init__(filename = filename, multispeciesdocket = multispeciesdocket)
        self.directory = directory
        self.path = self.directory + self.filename
        
