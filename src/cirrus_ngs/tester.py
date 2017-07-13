from dnaSeq import WGSPipelineManager
from cfnCluster import CFNClusterManager
from cfnCluster import ConnectionManager

master_ip_address = CFNClusterManager.create_cfn_cluster(cluster_name="mustafa")
ssh_client = ConnectionManager.connect_master(hostname=master_ip_address,username="ec2-user", private_key_file="/home/mustafa/interns_oregon_key.pem")

WGSPipelineManager.execute(ssh_client, "project_name", [""], "s3://ucsd-ccbb-interns/Mustafa/wgs_test/Sample_cDNA", [['Sample_cDNA93_R1.fq.gz', 'Sample_cDNA93_R2.fq.gz']], ['groupA'], "s3://ucsd-ccbb-interns/Mustafa/wgs_test/Sample_cDNA", "")
print("done")
