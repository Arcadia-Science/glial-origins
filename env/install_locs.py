##############################################################################
### Update these file paths to reflect your local directory structure ########

BASE_DIR = '/home/ec2-user/'
HOME_DIR = '/home/ec2-user/miniconda3/'
GFFREAD_LOC = HOME_DIR + 'envs/glial_origins/bin/gffread'
TRANSDECODER_LOC = HOME_DIR + 'pkgs/transdecoder-5.5.0-pl526_1/opt/transdecoder/'
TDLONGORF_LOC = HOME_DIR + 'pkgs/transdecoder-5.5.0-pl526_1/opt/transdecoder/TransDecoder.LongOrfs'
TDPREDICT_LOC = HOME_DIR + 'pkgs/transdecoder-5.5.0-pl526_1/opt/transdecoder/TransDecoder.Predict'
TDLONGESTORF_LOC = HOME_DIR + 'pkgs/transdecoder-5.5.0-pl526_1/opt/transdecoder/util/get_longest_ORF_per_transcript.pl'
ID_MAPPER_LOC = '../../utils/id_mapper.sh'
ID_MAPPER_GENENAME_LOC = '../../utils/id_mapper_genename.sh'
ORTHOFINDER_LOC = HOME_DIR + 'bin/orthofinder'

##############################################################################
##### The code below automatically detects the top of your git repository ####

import subprocess
GIT_HOME = subprocess.run(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).stdout.decode("utf-8").strip('\n')
GLOBAL_OUTPUT_DIRECTORY = GIT_HOME + '/output/'

##############################################################################
###### Change this address to point to your own AWS S3 bucket of choice ######

S3_BUCKET_ADDRESS = 's3://arcadia-reference-datasets/'

##############################################################################