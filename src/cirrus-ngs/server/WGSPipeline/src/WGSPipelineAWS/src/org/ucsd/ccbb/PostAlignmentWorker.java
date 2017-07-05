package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * PostAlignmentWorker processes realignment, variant calls and so on.
 * 
 * @author Guorong Xu
 */
public class PostAlignmentWorker {
	
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores the BAM file name?*/
	public String m_bamFile = "";
	/**a string stores the group name?*/
	public String m_group = "";
	/**a string stores the configuration file name?*/
	public String m_configFile = "";

	/**
	 * Constructs a PostAlignmentWorker instance with specified S3 upload URL, S3 
	 * download URL, BAM file, group and configuration file. 
	 * @param S3UploadURL		the S3 upload URL
	 * @param S3DownloadURL		the S3 download URL
	 * @param bamFile			the BAM file
	 * @param group				the group?
	 * @param configFile		the configuration file 
	 */
	public PostAlignmentWorker(String S3UploadURL, String S3DownloadURL, String bamFile, String group, String configFile) {
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_bamFile = bamFile;
		this.m_group = group;
		this.m_configFile = configFile;
	}

	/**
	 * Processes commands.
	 * First, download files from S3. Then perform post- alignment. Next, 
	 * print reads for variants calling. Then, uploads files to S3. Last, delete
	 * the workspace.
	 * If S3.delete.enable is true, then delete files on S3
	 * @throws Exception
	 */
	public void processCommand() throws Exception {		
		Date date= new Date();
		
		System.out.println("BAMWorker is processing " + m_bamFile);
		ConfigLoader loader = new ConfigLoader();
		loader.parsePathFile(m_configFile);
		
		String workspace = loader.get("work.space", "/scratch/workspace");	
		String logPath = loader.get("home.dirctory", "/shared/workspace/WGSPipeline") + "/results/logs";		
		File dir = new File(logPath);
		if (!dir.exists())
			dir.mkdirs();
		
		/*create an output file with the same name as the BAM file*/
		String outputFile = m_bamFile.substring(0, m_bamFile.indexOf(".sort.split.bam"));
        
		/*download files from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
	    S3DownloadExecutor s3d= new S3DownloadExecutor();
	    s3d.setShellFile(loader.get("s3download.shell", "/shared/workspace/WGSPipeline/scripts/s3download.sh"));
	    s3d.setFileName(m_bamFile);
	    s3d.setS3DownloadURL(m_S3DownloadURL);
	    s3d.setLogFile(logPath);
	    s3d.execute("", m_bamFile);
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
		
		//Performing post-alignment
		System.out.println(new Timestamp(date.getTime()) + ": Post-alignment is processing.");
		PostAlignmentExecutor pae = new PostAlignmentExecutor();
		pae.setShellFile(loader.get("post-alignment.shell", "/shared/workspace/WGSPipeline/scripts/postalignment.sh"));
		pae.setLogFile(logPath);
		pae.execute("", outputFile);
		System.out.println(new Timestamp(date.getTime()) + ": Post-alignment is done.");
		
		//Printing reads for variants calling
		System.out.println(new Timestamp(date.getTime()) + ": Print-reads is processing.");
		PrintReadsExecutor pre = new PrintReadsExecutor();
		pre.setShellFile(loader.get("post-alignment.print.reads.shell", "/shared/workspace/WGSPipeline/scripts/print_reads.sh"));
		pre.setFileName(outputFile);
		pre.setLogFile(logPath);
		pre.execute("", outputFile);
		System.out.println(new Timestamp(date.getTime()) + ": Print-reads is done.");
		
		/*upload files to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");
	    S3UploadExecutor s3u= new S3UploadExecutor();
	    s3u.setShellFile(loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setS3UploadURL(m_S3UploadURL);
	    s3u.setLogFile(logPath);
	    
	    s3u.setFileName(outputFile + ".final.bam");
	    s3u.execute("", outputFile + ".final.bam");
	    
	    s3u.setFileName(outputFile + ".final.bed");
	    s3u.execute("", outputFile + ".final.bed");
	    
	    s3u.setFileName(outputFile + ".final.bai");
	    s3u.execute("", outputFile + ".final.bai");
        System.out.println(new Timestamp(date.getTime()) + ": S3 upload is done.");
        
        /*delete workspace*/
	    System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is processing.");
	    DeleteWorkspaceExecutor deleter = new DeleteWorkspaceExecutor();
	    deleter.setInputFile(workspace + "/" + outputFile);
	    deleter.setShellFile(loader.get("delete.workspace.shell", "/shared/workspace/WGSPipeline/scripts/delete.sh"));
	    deleter.execute("", workspace + "/" + outputFile);
		System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is done.");
		
		/*delete files on S3*/
		if ("true".equalsIgnoreCase(loader.get("s3.delete.enabled", "false")))
		{
		    System.out.println(new Timestamp(date.getTime()) + ": S3 delete is processing.");
		    S3DeleteExecutor s3deleter = new S3DeleteExecutor();
		    s3deleter.setShellFile(loader.get("s3delete.shell", "/shared/workspace/WGSPipeline/scripts/s3delete.sh"));
		    s3deleter.setS3URL(m_S3DownloadURL);
		    s3deleter.setLogFile(logPath);
		    s3deleter.setFileName(m_bamFile);
		    s3deleter.execute("", m_bamFile);
			System.out.println(new Timestamp(date.getTime()) + ": S3 delete is done.");
		}     
	}
}