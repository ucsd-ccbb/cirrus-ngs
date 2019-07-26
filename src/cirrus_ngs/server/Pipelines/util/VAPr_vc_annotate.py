# Import core module
from VAPr import vapr_core
import os

# Start by specifying the project information
IN_PATH = "/path/to/vcf"
OUT_PATH = "/path/to/out"
ANNOVAR_PATH = "/path/to/annovar"
MONGODB = 'VariantDatabase'
COLLECTION = 'Cancer'

annotator = vapr_core.VaprAnnotator(input_dir=IN_PATH,
                                   output_dir=OUT_PATH,
                                   mongo_db_name=MONGODB,
                                   mongo_collection_name=COLLECTION,
                                   build_ver='hg19',
                                   vcfs_gzipped=False,
                                   annovar_install_path=ANNOVAR_PATH)
