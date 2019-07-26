from VAPr import vapr_core
from pymongo import MongoClient
import pandas as pd
import os
import sys
from VAPr.formatting import maf_formatter, create_whole_dataset_list

ALU_IN_PATH=sys.argv[1]
NONALU_IN_PATH=sys.argv[2]
OUT_PATH=sys.argv[3]
ANNOVAR_PATH=sys.argv[4]
MONGODB=sys.argv[5]
strand = sys.argv[6]

strand_annotation_df = pd.read_csv(strand, sep="\t")

# Alu
alu_annotator = vapr_core.VaprAnnotator(input_dir=ALU_IN_PATH,
					output_dir=OUT_PATH,
					mongo_db_name=MONGODB,
				        mongo_collection_name='alu',
					build_ver='hg19',
					vcfs_gzipped=True,
					annovar_install_path=ANNOVAR_PATH)

alu_annotator.annotate(num_processes=8)

client = MongoClient()
db = getattr(client, MONGODB)
alu_collection = getattr(db, 'alu')
alu_filtered = alu_collection.find({"$and": [{"dbsnp": {"$exists": False}},
				   {"1000g2015aug_all": {"$exists": False}},
				   {"$or":[{"$and": [{"ref": "A"}, {"alt": "G"}]}, {"$and": [{"ref": "T"},{"alt": "C"}]}]}]})

alu_filtered = list(alu_filtered)
alu_maf = maf_formatter(alu_filtered)
alu_maf["Novelty"] = "Alu"

# Nonalu
nonalu_annotator = vapr_core.VaprAnnotator(input_dir=NONALU_IN_PATH,
					   output_dir=OUT_PATH,
					   mongo_db_name=MONGODB,
					   mongo_collection_name='nonalu',
					   build_ver='hg19',
					   vcfs_gzipped=True,
					   annovar_install_path=ANNOVAR_PATH)

nonalu_annotator.annotate(num_processes=8)
db = getattr(client, MONGODB)
nonalu_collection = getattr(db, 'nonalu')
nonalu_filtered = nonalu_collection.find({"$and": [{"dbsnp": {"$exists": False}},
					{"1000g2015aug_all": {"$exists": False}},
				        {"$or":[{"$and": [{"ref": "A"}, {"alt": "G"}]}, {"$and": [{"ref": "T"},{"alt": "C"}]}]}]})

nonalu_filtered = list(nonalu_filtered)
nonalu_maf = maf_formatter(nonalu_filtered)
nonalu_maf["Novelty"] = "NonAlu"

# Combine
maf = pd.concat([alu_maf, nonalu_maf], sort=False)
maf["depth"] = maf.t_ref_count + maf.t_alt_count
maf["VAF"] = maf.t_alt_count/(maf.t_ref_count + maf.t_alt_count)
maf = maf[(maf['t_alt_count'] >= 2) & (maf['depth'] >= 5) & (maf['VAF'] >= 0.10) & (maf['VAF'] <= 0.95)]
maf = maf.merge(strand_annotation_df, on="Hugo_Symbol")
maf = maf[((maf.Reference_Allele == "A") & (maf.Strand == "+")) | ((maf.Reference_Allele == "T") & (maf.Strand == "-"))]

# Write out MAF
out_maf = os.path.join(OUT_PATH, MONGODB + '_vapr.maf')
maf.to_csv(out_maf, sep="\t", index=False)

