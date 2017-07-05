package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * GATKHaplotypeWorker applies HaplotypeCaller in GATK.
 * 
 * @author Guorong Xu
 */
public class GATKHaplotypeWorker {
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
	 * Constructs a GATKHaplotypeWorker instance with specified S3 upload URL,
	 * S3 download URL, BAM file, group name and configuration file. 
	 * @param S3UploadURL		the S3 upload URL
	 * @param S3DownloadURL		the S3 download URL
	 * @param bamFile			the BAM file
	 * @param group				the group
	 * @param configFile		the configuration file
	 */
	public GATKHaplotypeWorker(String S3UploadURL, String S3DownloadURL, String bamFile, String group, String configFile) {
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_bamFile = bamFile;
		this.m_group = group;
		this.m_configFile = configFile;
	}

	/**
	 * Processes commands.
	 * First, download files from S3. Then process GATKHaplotype. Next, 
	 * upload files to S3. Last, delete the workspace.
	 * If S3.delete.enable is true, then delete files on S3?
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
		String outputFile = m_bamFile.substring(0, m_bamFile.indexOf(".final"));
        
		/*download files from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
	    S3DownloadExecutor s3d= new S3DownloadExecutor();
	    s3d.setShellFile(loader.get("s3download.shell", "/shared/workspace/WGSPipeline/scripts/s3download.sh"));
	    s3d.setS3DownloadURL(m_S3DownloadURL);
	    s3d.setLogFile(logPath);
	    s3d.setFileName(m_bamFile + ".bam");
	    s3d.execute("", m_bamFile + ".bam");
	    
	    s3d.setFileName(m_bamFile + ".bai");
	    s3d.execute("", m_bamFile + ".bai");
	    
	    s3d.setFileName(m_bamFile + ".bed");
	    s3d.execute("", m_bamFile + ".bed");
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
		
        /*process GATKHaplotype*/
        System.out.println(new Timestamp(date.getTime()) + ": GATKHaplotype is processing.");
		GATKHaplotypeExecutor ghe = new GATKHaplotypeExecutor();
		ghe.setShellFile(loader.get("variants.calling.gatkhaplotype.shell", "/shared/workspace/WGSPipeline/scripts/gatkhaplotype.sh"));
		ghe.setFileName(outputFile);
		ghe.setGroup(m_group);
		ghe.setLogFile(logPath);
		ghe.execute("", outputFile);
		
		System.out.println(new Timestamp(date.getTime()) + ": GATKHaplotype is done.");
		
		/*upload files to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");
	    S3UploadExecutor s3u= new S3UploadExecutor();
	    s3u.setShellFile(loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setS3UploadURL(m_S3UploadURL);
	    s3u.setLogFile(logPath);
	    s3u.setFileName(outputFile + ".raw.vcf.gz");
	    s3u.execute("", outputFile + ".raw.vcf.gz");
	    
	    s3u.setFileName(outputFile + ".raw.vcf.gz.tbi");
	    s3u.execute("", outputFile + ".raw.vcf.gz.tbi");

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
		    s3deleter.setFileName(m_bamFile + ".bam");
		    s3deleter.execute("", m_bamFile + ".bam");
		    
		    s3deleter.setFileName(m_bamFile + ".bai");
		    s3deleter.execute("", m_bamFile + ".bai");
		    
		    s3deleter.setFileName(m_bamFile + ".bed");
		    s3deleter.execute("", m_bamFile + ".bed");
			System.out.println(new Timestamp(date.getTime()) + ": S3 delete is done.");
		}    
	}
}