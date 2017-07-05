import java.io.File;
import java.util.ArrayList;

public class SmallRNASeqExecutor {
	/**
	 * The Small RNASeq pipeline executor is to process the sequencing files and includes all modules.
	 * 
	 */
	public String m_shellScript = "";
	public ArrayList<String> m_fastq = new ArrayList<String>();
	public ArrayList<String> m_fastqFileName = new ArrayList<String>();
	
	public String m_suffix = "";
		
	public static void main(String[] args) throws Exception 
	{
		if (args.length < 1)
		{
			System.out.println("Usage: java -jar SmallRNASeqPipeline.jar <yaml_file>");
			System.out.println("For example: java -jar SmallRNASeqPipeline.jar /path/to/yaml_file");
			
			return;
		}
		
		String rootPath = "/shared/workspace/SmallRNASeqPipeline/";
		String yamlFile = args[0];
		/*create a new YamlParser instance*/
		YamlParser parser = YamlParser.getInstance();
		/*parse the yaml file*/
		parser.parse(yamlFile);
		
		String analysis = parser.getAnalysis();	
		SmallRNASeqExecutor executor = new SmallRNASeqExecutor();
		
		int nodeNum = new PBSTracker().checkNodeNum();
		
		System.out.println();
		System.out.println("Cluster nodes number is: " + nodeNum);
		System.out.println();
		
		//The executor processes the pipeline.
		executor.processFiles( rootPath,  yamlFile,  analysis);
	}
	 
    /**
     * The method is to process sequencing files.
     * @param loader
     * @param yamlFile
     * @param inputFolder
     * @param suffix
     * @throws Exception
     */
	public void processFiles(String rootPath, String yamlFile, String analysis) throws Exception
	{
		/*create a new ConfigLoader instance*/
		ConfigLoader loader = new ConfigLoader();
		/*process the system.conf file*/
		loader.parsePathFile(rootPath + "system.conf");
		
		YamlParser parser = YamlParser.getInstance();
		String projectInfo = parser.m_projectName;
		String S3UploadURL = parser.getS3UploadURL();
		String analysisSteps = parser.getAnalysis();
		ArrayList<String[]> fastqList = parser.getFastqList();
		
		System.out.println("Total number of FASTQ files: " + fastqList.size() + ".");
		System.out.println();
		
		String dataFolder = "/shared/workspace/data_archive/SmallRNASeq/" + projectInfo + "/";
	
		// To call S3DownloadExecutor to download fastq files.
		S3DownloadExecutor s3download = new S3DownloadExecutor();
					
	    for ( int index = 0; index < fastqList.size(); index++ ) 
	    {
			String fileName = fastqList.get(index)[0];
			String S3DownloadURL = fastqList.get(index)[1];
					
			File dir = new File(dataFolder);
			if (!dir.exists())
				dir.mkdirs();
			
			System.out.println("System is downloading: " + fileName);
			s3download.setDataFolder(dataFolder);
			s3download.setS3DownloadURL(S3DownloadURL);
			s3download.run(loader.m_configsList.get("s3.download.shell"), fileName);
		}	
		new PBSTracker().trackPBSQueue(Integer.valueOf(loader.m_configsList.get("trackTime")), "download");
		
		if (analysisSteps.indexOf("fastqc") > -1)
		{
			// To call fastQC module to run fastQC.
			FastQCExecutor fe = new FastQCExecutor();
						
		    for ( int index = 0; index < fastqList.size(); index++ ) 
		    {
				String fileName = fastqList.get(index)[0];
						
				File dir = new File(dataFolder);
				if (!dir.exists())
					dir.mkdirs();
				
				String fastqFileName = dataFolder + fileName;	
				System.out.println("System is processing: " + fastqFileName);
				fe.run(loader.m_configsList.get("fastqc.shell"), fastqFileName);
			}	
			new PBSTracker().trackPBSQueue(Integer.valueOf(loader.m_configsList.get("trackTime")), "fastQC");
		}
		
		if (analysisSteps.indexOf("bowtie-alignment") > -1)
		{
			// To run alignment using Bowtie.
			AlignmentExecutor ne = new AlignmentExecutor();
		    for ( int index = 0; index < fastqList.size(); index++ ) 
		    {
				String fileName = fastqList.get(index)[0];
				
				File dir = new File(dataFolder);
				if (!dir.exists())
					dir.mkdirs();
				
				String fastqFileName = dataFolder + fileName;		
				ne.runBowtie(loader.m_configsList.get("alignment.bowtie.se.shell"), dataFolder, fastqFileName);		
			}	
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.m_configsList.get("trackTime")), "alignment");
		    
			//To merge the alignment log files.
			MergeLogsExecutor mle = new MergeLogsExecutor();
			mle.run(loader.m_configsList.get("mergeLogs.shell"), dataFolder);
			new PBSTracker().trackPBSQueue(Integer.valueOf(loader.m_configsList.get("trackTime")), "mergeLogs");
		}	
	
		if (analysisSteps.indexOf("counting") > -1)
		{
			MiRNAParser mirParser = new MiRNAParser();		
			mirParser.parseMiRBase(dataFolder, dataFolder + "output/", ".sam", loader.m_configsList.get("mirbase.annotation"), yamlFile, "true");
		}
		
		System.out.println("");
		// To call S3UploadExecutor to upload files.
		S3UploadExecutor s3upload = new S3UploadExecutor();		
		
		s3upload.setDataFolder(dataFolder);
		s3upload.setS3UploadURL(S3UploadURL);
		s3upload.run(loader.m_configsList.get("s3.upload.shell"), dataFolder);
	
		new PBSTracker().trackPBSQueue(Integer.valueOf(loader.m_configsList.get("trackTime")), "upload");
		
		System.out.println("");
		System.out.println("");
		System.out.println("======================================================");
		System.out.println("The processing of the project \"" + projectInfo + "\" is done!");
		System.out.println("======================================================");
}
	
	/**
	 * The method is to check if the output file output correctly,
	 * otherwise we need to rerun it again.
	 */
	public void checkOutputFiles()
	{
		ArrayList<String> fastq = new ArrayList<String>();
		ArrayList<String> fastqFileName = new ArrayList<String>();
		
		for (int i = 0; i < m_fastq.size(); i++)
		{
			String fastqFile = m_fastq.get(i);
			String outputSAMFile = fastqFile.substring(0, fastqFileName.indexOf("1."));
			if (!(new File(outputSAMFile + ".sam").exists() && new File(outputSAMFile + ".bam").exists()
					&& new File(outputSAMFile + ".bam.bai").exists()))
			{
				fastq.add(m_fastq.get(i));
				fastqFileName.add(m_fastqFileName.get(i));
			}
		}
		
		m_fastq = fastq;
		m_fastqFileName = fastqFileName;
	}
}
