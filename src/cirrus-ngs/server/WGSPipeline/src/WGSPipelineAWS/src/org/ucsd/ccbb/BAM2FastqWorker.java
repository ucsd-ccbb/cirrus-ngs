package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * BAM2FastqWorker converts BAM format to Fastq format.
 * 
 * @author Guorong Xu
 */
public class BAM2FastqWorker {
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores the name of the BAM File*/
	public String m_bamFile = "";
	/**a string stores the name of the configuration file*/
	public String m_configFile = "";
	
	/**
	 * Construct a new BAM2FastqWorker instance with specified S3 upload URL, 
	 * S3 download URL, BAM file name and configuration file name.
	 * @param S3UploadURL	S3 upload URL
	 * @param S3DownloadURL	S3 download URL
	 * @param bamFile		BAM file name
	 * @param configFile	configuration file name
	 */
	public BAM2FastqWorker(String S3UploadURL, String S3DownloadURL, String bamFile, String configFile) {
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_bamFile = bamFile;
		this.m_configFile = configFile;
	}

	/**
	 * The function processes commands.
	 * First, download files from S3. Then convert BAM format to fastq format.
	 * Next, compress fastq format files. Last, upload files to S3.
	 * @throws Exception
	 */
	public void processCommand() throws Exception {		
		Date date= new Date();
		
		System.out.println(new Timestamp(date.getTime()) + ": BAM2FastqWorker is processing " + m_bamFile);
		ConfigLoader loader = new ConfigLoader();
		loader.parsePathFile(m_configFile);
		
		String workspace = loader.get("work.space", "/scratch/workspace");		
		String logPath = loader.get("home.dirctory", "/shared/workspace/WGSPipeline") + "/results/logs";	
		System.out.println("logPath: " + logPath);
		File dir = new File(logPath);
		if (!dir.exists())
			dir.mkdirs();
		
		/*download the file from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
	    S3DownloadExecutor s3d= new S3DownloadExecutor();
	    s3d.setShellFile(loader.get("s3download.shell", "/shared/workspace/WGSPipeline/scripts/s3download.sh"));
	    s3d.setFileName(m_bamFile);
	    s3d.setS3DownloadURL(m_S3DownloadURL);
	    s3d.setLogFile(logPath);
	    s3d.execute("",  m_bamFile);
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
        
        /*convert BAM format to fastq format*/
	    System.out.println(new Timestamp(date.getTime()) + ": BAM to Fastq is processing.");
		BAM2FastqExecutor b2f = new BAM2FastqExecutor();
		b2f.setShellFile(loader.get("bam.to.fastq.shell", "/shared/workspace/WGSPipeline/scripts/bam_to_fastq.sh"));
		b2f.setLogFile(logPath);
		b2f.execute("", m_bamFile.substring(0, m_bamFile.indexOf(".")));
        System.out.println(new Timestamp(date.getTime()) + ": BAM to Fastq is done.");
             
        /*compress fastq file*/
        System.out.println(new Timestamp(date.getTime()) + ": Compress-Fastq is processing.");	
        /* creates a thread pool that reuses 2 threads operating off a shared 
         * unbounded queue*/
		ExecutorService es = Executors.newFixedThreadPool(2);		
		for (int i = 1; i < 3; i ++)
		{			
			String fastqFile = m_bamFile.substring(0, m_bamFile.indexOf(".")) + "_R" + i + ".fq";
			Runnable worker = new ZipFastqWorker(loader, fastqFile, logPath);
			/* Executes the given command*/
			es.execute(worker);
		}
		/* Shut down es.
		 * Allow previously submitted tasks to execute before terminating*/
		es.shutdown();
		while (!es.isTerminated()) {}
		System.out.println(new Timestamp(date.getTime()) + ": Compress-Fastq is done.");
		
		/*upload fastq format files to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");	
		es = Executors.newFixedThreadPool(2);		
		for (int i = 1; i < 3; i ++)
		{			
			String fastqFile = m_bamFile.substring(0, m_bamFile.indexOf(".")) + "_R" + i + ".fq.gz";
			Runnable worker = new S3UploadWorker(loader, m_S3UploadURL, fastqFile, logPath);
			es.execute(worker);
		}
		
		es.shutdown();
		while (!es.isTerminated()) {}	    
        System.out.println(new Timestamp(date.getTime()) + ": S3 upload is done.");
	}
}