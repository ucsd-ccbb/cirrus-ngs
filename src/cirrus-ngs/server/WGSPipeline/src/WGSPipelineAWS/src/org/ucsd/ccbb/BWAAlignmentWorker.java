package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * This class processes all alignments, sorting and splitting.
 *
 * @author Guorong Xu
 */
public class BWAAlignmentWorker {
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores the name of the first fastqFile*/
	public String m_fastqFile1 = "";
	/**a string stores the name of the second fastqFile*/
	public String m_fastqFile2 = "";
	/**a string stores the name of the configuration file*/
	public String m_configFile = "";
	/**a string stores the description*/
	public String m_description = "";
	/**a string stores the date*/
	public String m_date = "";
	
	/**
	 * Construct a new BWAAlignmentWorker instance with specified S3 upload 
	 * URL, S3 download URL, the first fastq file name, the second fastq file
	 * name, description, date and configuration file name.
	 * @param S3UploadURL		S3 upload URL
	 * @param S3DownloadURL		S3 download URL
	 * @param fastqFile1		the first fastq file name
	 * @param fastqFile2		the second fastq file name
	 * @param description		description of the project
	 * @param date				date of the project
	 * @param configFile		configuration file name
	 */
	public BWAAlignmentWorker(String S3UploadURL, String S3DownloadURL, String fastqFile1, String fastqFile2, 
			String description, String date, String configFile) 
	{
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_fastqFile1 = fastqFile1;
		this.m_fastqFile2 = fastqFile2;
		this.m_description = description;
		this.m_date = date;
		this.m_configFile = configFile;
	}

	/**
	 * The function processes commands.
	 * First, download files from S3. Then perform alignment with BWA. Next, 
	 * sort all BAM files. Then, split a BAM file into small BAM files. Next,
	 * list all splitted BAM files. Last, delete the workspace. 
	 * @throws Exception
	 */
	public void processCommand() throws Exception {		
		Date date= new Date();
		
		System.out.println(new Timestamp(date.getTime()) + ": System is processing " + m_fastqFile1 + "," + m_fastqFile2);
		ConfigLoader loader = new ConfigLoader();
		loader.parsePathFile(m_configFile);
		
		String workspace = loader.get("work.space", "/scratch/workspace");	
		String logPath = loader.get("home.dirctory", "/shared/workspace/WGSPipeline") + "/results/logs";		
		File dir = new File(logPath);
		if (!dir.exists())
			dir.mkdirs();
		
		/*create an output file with the same name as fastqFiles*/
		String outputFile = m_fastqFile1.substring(0, m_fastqFile1.indexOf("_R"));
        
		/*download files from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
	    ExecutorService es = Executors.newFixedThreadPool(2);	
        for (int i = 1; i < 3; i++)
        {
 			Runnable worker = new S3DownloadWorker( loader, m_S3DownloadURL, outputFile + "_R" + i + ".fq.gz", logPath);
			es.execute(worker);
	    }		
		es.shutdown();
		while (!es.isTerminated()) {}
        
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
        
		// Performing alignment with BWA
		System.out.println(new Timestamp(date.getTime()) + ": Alignment is processing.");
		BWAExecutor bwa = new BWAExecutor();
		bwa.setShellFile(loader.get("alignment.bwa.pe.shell", "/shared/workspace/WGSPipeline/scripts/bwa.sh"));			
		bwa.setDate(m_date);
		bwa.setDescription(m_description);
		bwa.setFastqFiles(m_fastqFile1, m_fastqFile2);
		bwa.setLogFile(logPath);
	 	bwa.execute("", outputFile);
        System.out.println(new Timestamp(date.getTime()) + ": Alignment processing is done.");
        
        //Sorting BAM file
	    System.out.println(new Timestamp(date.getTime()) + ": Sorting BAM is processing.");
		SortBAMExecutor sorter = new SortBAMExecutor();
		sorter.setInputFile(outputFile);
		sorter.setLogFile(logPath);
		sorter.setShellFile(loader.get("post-alignment.sort.bam.shell", "/shared/workspace/WGSPipeline/scripts/sort_bam.sh"));
		sorter.execute("", outputFile);
		System.out.println(new Timestamp(date.getTime()) + ": Sorting BAM is done.");
       
		//Splitting BAM file
        System.out.println(new Timestamp(date.getTime()) + ": Splitting BAM is processing.");
		es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));	
        ChromSplitter splitter = new ChromSplitter();
        ArrayList<String[]> regions = splitter.getRegions(loader);
        for (int i = 0; i < regions.size(); i++)
        {
        	String[] region = regions.get(i);
			//System.out.println("Splitting sam/bam file: " + chrom + ":" + start + "-" + end);			
			Runnable worker = new SplitWorker( loader, outputFile, region[0], Integer.valueOf(region[1]), Integer.valueOf(region[2]), logPath);
			es.execute(worker);
        }
		
		es.shutdown();
		while (!es.isTerminated()) {}
        System.out.println(new Timestamp(date.getTime()) + ": Splitting BAM is done.");
        
		//List all split BAM files
		PipelineUtil util = new PipelineUtil();
		util.listAllFiles(workspace, ".sort.split.bam");
		ArrayList<String> splitBAMs = util.getFileNameList();
		
		/*Upload all sorted and splitted BAM files to S3*/
		System.out.println(new Timestamp(date.getTime()) + ": Uploading all sorted and split bam files to S3.");
		es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));		
		for (int i = 0; i < splitBAMs.size(); i ++)
		{			
			Runnable worker = new S3UploadWorker(loader, m_S3UploadURL + "/" + outputFile, splitBAMs.get(i), logPath);
			es.execute(worker);
		}
		
		es.shutdown();
		while (!es.isTerminated()) {}
		
		//Deleting workspace
	    System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is processing.");
	    DeleteWorkspaceExecutor deleter = new DeleteWorkspaceExecutor();
	    deleter.setInputFile(workspace);
	    deleter.setShellFile(loader.get("delete.workspace.shell", "/shared/workspace/WGSPipeline/scripts/delete.sh"));
	    deleter.execute("", workspace);
		System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is done.");       
	}
}