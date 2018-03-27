# cirrus-ngs

cloud-optimized primary analysis pipelines for rna-seq, mirna-seq, chip-seq, and variant calling in whole-genome/whole-exome dna-seq.

## introduction

bionformatic analysis of large-scale next-generation sequencing (ngs) data requires significant compute resources. while cloud computing makes such processing power available on-demand, the administration of dynamic compute clusters is a daunting task for most working biologists.  to address this pain point, the center for computational biology and bioinformatics at the university of california, san diego, has developed cirrus-ngs, a turn-key solution for common ngs analyses using amazon web services (aws).
	
cirrus-ngs provides primary analysis pipelines for rna-seq, mirna-seq, chip-seq, and whole-genome/whole-exome sequencing data.  cirrus users need not have any bioinformatics tools for these pipelines installed on their local machine, as all computation is performed using aws clusters and all results are uploaded to aws's s3 remote storage.  the clusters are created dynamically through cirrus-ngs, based on a custom machine image with all the necessary software pre-installed.  users manage clusters and run pipelines from within a web browser, using flexible but light-weight jupyter notebooks.

## installation

### getting an aws account

because cirrus-ngs employs aws for computation and storage, users must have an active aws account.  if you do not have such an account, visit [https://docs.aws.amazon.com/awsec2/latest/userguide/get-set-up-for-amazon-ec2.html](https://docs.aws.amazon.com/awsec2/latest/userguide/get-set-up-for-amazon-ec2.html) for guidance on how to create one, and execute the steps it describes. 


### installing `conda`

while cirrus-ngs itself is provied through github, its supporting libraries are best installed through conda, a cross-platform package manager.  if you do not have conda installed already, download the appropriate python-3.6-based installer from the links below:

* linux: [https://repo.continuum.io/miniconda/miniconda3-latest-linux-x86_64.sh](https://repo.continuum.io/miniconda/miniconda3-latest-linux-x86_64.sh)
* osx: [https://repo.continuum.io/miniconda/miniconda3-latest-macosx-x86_64.sh](https://repo.continuum.io/miniconda/miniconda3-latest-macosx-x86_64.sh)

then run the installer with the following script, replacing `<miniconda_py3.sh>` with the name of the install script downloaded above:

	bash <miniconda_py3.sh> -b -p $home/miniconda3
	echo "export path=\"$home/miniconda3/bin:\$path\"" >>$home/.bashrc
	source $home/.bashrc	

### installing cirrus-ngs

cirrus-ngs is available for linux-64 or osx-64 platforms, and requires python 3.5 or higher.

in the directory in which you would like to install cirrus-ngs, run the following commands:

	conda install paramiko pyyaml git jupyter notebook
	pip install scp cfncluster
	git clone https://github.com/ucsd-ccbb/cirrus-ngs.git

these commands install the necessary supporting libraries and create a new directory called `cirrus-ngs` that holds the cirrus software.

### starting the notebook server

since the cirrus-ngs interface is provided through jupyter notebooks, the first step of running cirrus is to start a local jupyter notebook server. from the `cirrus-ngs` directory, run:

	cd notebooks/cirrus-ngs
	jupyter notebook

you will receive a message like

	[i 10:48:39.779 notebookapp] serving notebooks from local directory: /users/<yourname>/cirrus-ngs/notebooks
	[i 10:48:39.779 notebookapp] 0 active kernels 
	[i 10:48:39.779 notebookapp] the jupyter notebook is running at: http://localhost:8889/?token=8dcb989e6852e3dfd679307470e17696c32771432c881573
	[i 10:48:39.779 notebookapp] use control-c to stop this server and shut down all kernels (twice to skip confirmation).
	
**note the url at which the server is running**.  it is usually http://localhost:8888, but if your system already has something running on port 8888, it may default so something else (e.g., in the example above, it is http://localhost:8889).

**do not shut down the jupyter notebook server** or close the terminal window in which it is running until you are done using cirrus-ngs, as it must be running in order for the cirrus-ngs notebooks to function.  when you are done using cirrus and ready to shut down the notebook server, simply type `control-c` in the terminal window where it is running and follow the prompts.

## running cirrus-ngs

### general overview

a user begins by creating a cluster.  after this, sh/e decides upon the pipeline to run, which determines the exact inputs required, and then creates (in a spreadsheet program or text editor) a design file that specifies the details of the input sequencing data.  s/he starts the jupyter notebook for the chosen pipeline and inputs the path to the design file and other required inputs (such as aws credentials, etc), then executes the notebook code to start the pipeline processing.  

the notebook creates a yaml file summarizing all of the user input and transfers that file to the cluster, where it directs cluster-native code to sequentially execute the analysis steps specified by the user in a distributed fashion. upon completion of each step, the output is uploaded to the user's specified s3 output bucket, from which it can be accessed at any point. during processing, the user can check the status of the pipeline run and manage the cluster status from within the notebook.

**a note about terminology**: throughout, the term "pipeline" is used to describe functionality for analyzing a general type of sequencing (such as rna-seq or whole-genome sequencing).  each pipeline can have multiple, more granular "workflows" that support different approaches to processing the same type of sequencing data; for instance, the rna-seq pipeline contains four different workflows (star_gatk, star_htseq, star_rsem, and kallisto).

### gathering required inputs

before beginning work with cirrus-ngs, gather the following information:

1. the aws access key id
	* see [https://docs.aws.amazon.com/general/latest/gr/managing-aws-access-keys.html](https://docs.aws.amazon.com/general/latest/gr/managing-aws-access-keys.html) to find your existing access key or create a new one)
2. the aws secret access key (see above link for more information)
3. a valid key pair (.pem) file for the account
	* if you encounter a `permission denied (publickey) error`, note that the permissions on your key (.pem) file **must** be set so that it is not public, e.g. by running `chmod 0400 <mypemfilename>.pem` in the directory where the .pem file is located.
	* see [https://docs.aws.amazon.com/awsec2/latest/userguide/ec2-key-pairs.html#having-ec2-create-your-key-pair](https://docs.aws.amazon.com/awsec2/latest/userguide/ec2-key-pairs.html#having-ec2-create-your-key-pair) for further details.
<!--todo: add region-->
4. the vpc (virtual private cloud) id for the account
	* if you do not know the vpc id, you can find it by going part-way through the process of launching an aws compute instance.  follow the instructions at [https://docs.aws.amazon.com/awsec2/latest/userguide/launching-instance.html](https://docs.aws.amazon.com/awsec2/latest/userguide/launching-instance.html) up until step 6 (it does not matter what values you select during these steps, as you won't need to actually launch the instance).  step 6 brings you to the "configure instance details" screen, on which you can find your vpc id, as well as your subnet ids (see next item):
	
	![configure instance details screen](docs/configure_instance_details.png)
	
5. the subnet id for the account for the region in which you wish to create the cluster (see above for how to find this value)
4. the url of an existing aws s3 remote storage bucket where you have placed your sequencing data.
	* visit [https://docs.aws.amazon.com/amazons3/latest/gsg/getstartedwiths3.html](https://docs.aws.amazon.com/amazons3/latest/gsg/getstartedwiths3.html) and execute the "sign up for amazon s3" and "create a bucket" steps (click on the boxes in the workflow diagram for detailed instructions) if you do not already have a bucket set up.
	* to place your sequencing files into the bucket, follow the directions for the "add an object to a bucket" step, or use an s3-enabled transfer client like [cyberduck](https://cyberduck.io/).
5. (optional) the url of a second existing aws s3 remote storage bucket for cirrus-ngs outputs.

### creating a cluster

run the following steps to create a cluster.  you may then use that cluster for all cirrus-ngs work going forward without having to rerun these steps (although you may also rerun them to create additional clusters, if desired).

1. in your browser, visit the address at which the jupyter notebook server is running
	* this is usually http://localhost:8888, but see the [starting the notebook server](#starting-the-notebook-server) section for additional details.
	* you will see a list of all the cirrus-ngs notebooks:

	![cirrus-ngs notebooks list](docs/notebooks_list.png) 
2. click on the `basiccfnclustersetup.ipynb` notebook to start it.
3. fill in the parameters in the first two cells in the notebook, using the values identified above in the [gathering required inputs](#gathering-required-inputs) section where relevant. <!--todo:expand-->
4. run the notebook.<!--todo:expand-->


### running a pipeline

1. choose the pipeline you wish to run.
	* cirrus-ngs currently offers pipelines for rna-seq, mirna-seq, chip-seq, and variant calling in whole-genome/whole-exome dna-seq.
2. create a design file for your data.
	* the design file is a tab-separated text file with either two or three columns (depending on the workflow chosen) that specifies the names of the sequence files to process and the necessary metadata describing them.
	* see the [building a design file](#building-a-design-file) section below for full specifications of the design file format.
3. in your browser, visit the address at which the jupyter notebook server is running
	* this is usually http://localhost:8888, but see the [starting the notebook server](#starting-the-notebook-server) section for additional details.
	* you will see a list of all the cirrus-ngs notebooks:

	![cirrus-ngs notebooks list](docs/notebooks_list.png)
4. click on the appropriate notebook for your chosen pipeline to start it.
5. fill in the parameters in the first two cells in the notebook, using the values identified above in the [gathering required inputs](#gathering-required-inputs) section where relevant.<!--todo:expand-->
6. run the notebook.<!--todo:expand-->
7. after the pipeline processing is complete, download the output files from your s3 bucket for further analysis locally.
	* to download the outputs from your bucket, follow the directions for the "view an object" step at [https://docs.aws.amazon.com/amazons3/latest/gsg/openinganobject.html](https://docs.aws.amazon.com/amazons3/latest/gsg/openinganobject.html) or use an s3-enabled transfer client like [cyberduck](https://cyberduck.io/).

## building a design file
the design file is the primary user input to cirrus-ngs. it is a tab-separated text file that specifies the names of the sequence files to process and the necessary metadata describing them.  for the rna-seq and mirna-seq pipelines, it contains two columns, while for the chip-seq and whole-exome/whole-genome sequencing pipelines, it contains three columns (of which the first two are the same as in the two-column format).  the design file has no header line.

### two-column format
in the two-column format, the **first column** is the filename of the sample (with extensions: e.g. fastq.gz), while the **second column** is the name of the group associated with that sample.  (group names are up to the user, but they are generally set to experimental conditions. for example, all control samples can be given a group named "control".)

for example, a two-column design file for two single-end-sequenced samples named `mysample1` and `mysample2` might look like:

```
	mysample1.fastq.gz		groupa
	mysample2.fastq.gz		groupb
```

if the sequencing data is paired-end, the first column includes the name of the forward sequencing file, followed by a comma, followed by the name of the reverse sequencing file (note that there must **not** be any spaces between these two file names--only a comma!)  an example two-column design file for two paired-end-sequenced samples named `mysample1` and `mysample2` might look like:

```
	mysample1_fwd.fastq.gz,mysample1_rev.fastq.gz		groupa
	mysample2_fwd.fastq.gz,mysample2_rev.fastq.gz		groupb
```

### three-column format
the three-column format has the same first two columns as the two-column format, but adds a third column containing a pipeline-specific identifier that differentiates sample types for the chip-seq and whole-genome/whole-exome sequencing pipelines.

* for the chip-seq pipeline, each sample is identified as either chip or input.
* for the whole-genome/whole-exome sequencing pipeline, each sample is identified as either tumor or normal.

the following constraints apply to three-column design files:
	
* the sample type identifiers are **case-sensitive**.
	* for example, a three-column design file for two paired-end-sequenced chip-seq samples named `mysample1` and `mysample2` might look like 
	
	```
	mysample1_fwd.fastq.gz,mysample1_rev.fastq.gz		groupa		chip
	mysample2_fwd.fastq.gz,mysample2_rev.fastq.gz		groupb		input
	```

* each three-column design file must use **either** chip and input sample type identifiers **or** normal and tumor sample type identifiers, but not both.
	* a design file like the example below would thus be **invalid**:

	```
	mysample1_fwd.fastq.gz,mysample1_rev.fastq.gz		badexample1		chip
	mysample2_fwd.fastq.gz,mysample2_rev.fastq.gz		badexample1		tumor
	```
* each group must have exactly one sample designed with each of the two relevant sample type identifiers.	* a design file like the example below would thus be **invalid**:

	```
	mysample1_fwd.fastq.gz,mysample1_rev.fastq.gz		badexample2		normal
	mysample2_fwd.fastq.gz,mysample2_rev.fastq.gz		badexample2		normal
	```	
* if two samples are designated as forming a chip/input or normal/tumor pair, they must both have the same group.
	* a design file like the example below would thus be **invalid**:

	```
	mysample1_fwd.fastq.gz,mysample1_rev.fastq.gz		badexample3		normal
	mysample2_fwd.fastq.gz,mysample2_rev.fastq.gz		badexample4		tumor
	```
