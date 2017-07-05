package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * MergeWorker merges splitted final BAM files and VCF files.
 * 
 * @author Guorong Xu
 */
public class MergeWorker {
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores the file name?*/
	public String m_fileName = "";
	/**a string stores the configuration file name?*/
	public String m_configFile = "";
	
	/**
	 * Constructs a MergeWorker instance with specified S3 upload URL, S3 
	 * download URL, file name and configuration file. 
	 * @param S3UploadURL		the S3 upload URL
	 * @param S3DownloadURL		the S3 download URL
	 * @param fileName			the file name
	 * @param configFile		the configuration file
	 */
	public MergeWorker(String S3UploadURL, String S3DownloadURL, String fileName, String configFile) {
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_fileName = fileName;
		this.m_configFile = configFile;
	}

	/**
	 * Merges small BAM files and VCF files.
	 * If S3.delete.enable is true, then delete files on S3
	 * @throws Exception
	 */
	public void processCommand() throws Exception {		
		Date date= new Date();
		
		System.out.println("MergeWorker is processing " + m_fileName);
		ConfigLoader loader = new ConfigLoader();
		loader.parsePathFile(m_configFile);
		
		String workspace = loader.get("work.space", "/scratch/workspace");	
		String logPath = loader.get("home.dirctory", "/shared/workspace/WGSPipeline") + "/results/logs";		
		File dir = new File(logPath);
		if (!dir.exists())
			dir.mkdirs();
		      
		/*download small BAM files from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
		ExecutorService es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));	
        ChromSplitter splitter = new ChromSplitter();
        ArrayList<String[]> regions = splitter.getRegions(loader);
        ArrayList<String> bamFileList = new ArrayList<String>();
        for (int i = 0; i < regions.size(); i++)
        {
        	String[] region = regions.get(i);
        	bamFileList.add(m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".final.bam");
			Runnable downloadWorker = new S3DownloadWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".final.bam", logPath);
			es.execute(downloadWorker);
	    }
		
		es.shutdown();
		while (!es.isTerminated()) {}
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
		
        /* Merge small BAM files*/
	    System.out.println(new Timestamp(date.getTime()) + ": Merge-BAM-files is processing.");
	    MergeBAMExecutor mbe= new MergeBAMExecutor();
	    mbe.setShellFile(loader.get("merge.bam.shell", "/shared/workspace/WGSPipeline/scripts/merge_bam.sh"));
	    mbe.setInputFileNames(bamFileList);
	    mbe.setOutputFileName(m_fileName);
	    mbe.setLogFile(logPath);
	    mbe.execute("", m_fileName);
        System.out.println(new Timestamp(date.getTime()) + ": Merge-BAM-files is done.");
		
        /*upload merged BAM file to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");
	    S3UploadExecutor s3u= new S3UploadExecutor();
	    s3u.setShellFile(loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setS3UploadURL(m_S3DownloadURL);
	    s3u.setLogFile(logPath);
	    s3u.setFileName(m_fileName + ".final.bam");
	    s3u.execute("", m_fileName + ".final.bam");
        System.out.println(new Timestamp(date.getTime()) + ": S3 upload is done.");
        
        /*delete workspace*/
	    System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is processing.");
	    DeleteWorkspaceExecutor deleter = new DeleteWorkspaceExecutor();
	    deleter.setInputFile(workspace + "/" + m_fileName);
	    deleter.setShellFile(loader.get("delete.workspace.shell", "/shared/workspace/WGSPipeline/scripts/delete.sh"));
	    deleter.execute("", workspace + "/" + m_fileName);
		System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is done.");
      
		/*download small VCF files from S3*/
		System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
		es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));	
        for (int i = 0; i < regions.size(); i++)
        {
        	String[] region = regions.get(i);
			Runnable downloadWorker = new S3DownloadWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".raw.vcf.gz", logPath);
			es.execute(downloadWorker);
	    }
		
		es.shutdown();
		while (!es.isTerminated()) {}
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
        
        /* Merge small VCF files by regions*/
	    System.out.println(new Timestamp(date.getTime()) + ": MergeVCF by regions is processing.");
		es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));		
		ArrayList<String> chroms = splitter.getChoms();
		for (int i = 0; i < chroms.size(); i++)
		{
			ArrayList<String> vcfFileList = new ArrayList<String>();
			for (int j = 0; j < regions.size(); j++)
	        {
	        	String[] region = regions.get(j);
	        	if (region[0].equalsIgnoreCase(chroms.get(i)))
	        		vcfFileList.add(m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".raw.vcf.gz");	    
	        }
			
			Runnable mergeVCFWorker = new MergeVCFWorker( loader, vcfFileList, m_fileName + "." + chroms.get(i),  "true", logPath);
			es.execute(mergeVCFWorker);	
		}
		es.shutdown();
		while (!es.isTerminated()) {}
		System.out.println(new Timestamp(date.getTime()) + ": MergeVCF by regions is done.");
		
		/*Merge VCF files by names of chromosomes*/
		System.out.println(new Timestamp(date.getTime()) + ": MergeVCF by chroms is processing.");
		ArrayList<String> vcfFileList = new ArrayList<String>();
		for (int i = 0; i < chroms.size(); i++)
			vcfFileList.add(m_fileName + "." + chroms.get(i) + ".raw.vcf.gz");		    
		
		es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));	
		Runnable mergeVCFWorker = new MergeVCFWorker( loader, vcfFileList, m_fileName,  "false", logPath);
		es.execute(mergeVCFWorker);
		es.shutdown();
		while (!es.isTerminated()) {}
        System.out.println(new Timestamp(date.getTime()) + ": MergeVCF by chroms is done.");
		
        /*upload merged VCF  vfeededc files to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");	    
	    s3u= new S3UploadExecutor();
	    s3u.setShellFile(loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setS3UploadURL(m_S3UploadURL);
	    s3u.setLogFile(logPath);
	    s3u.setFileName(m_fileName + ".final.vcf.gz");
	    s3u.execute("", m_fileName + ".final.vcf.gz");
	    s3u.setFileName(m_fileName + ".final.vcf.gz.tbi");
	    s3u.execute("", m_fileName + ".final.vcf.gz.tbi");
        System.out.println(new Timestamp(date.getTime()) + ": S3 upload is done.");
        
        /*delete workspace*/
	    System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is processing.");
	    deleter = new DeleteWorkspaceExecutor();
	    deleter.setInputFile(workspace + "/" + m_fileName);
	    deleter.setShellFile(loader.get("delete.workspace.shell", "/shared/workspace/WGSPipeline/scripts/delete.sh"));
	    deleter.execute("", workspace + "/" + m_fileName);
		System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is done.");   
		
		/*delete files in S3*/
		if ("true".equalsIgnoreCase(loader.get("s3.delete.enabled", "false")))
		{
		    System.out.println(new Timestamp(date.getTime()) + ": S3 delete is processing.");
			es = Executors.newFixedThreadPool(Integer.valueOf(loader.get("slots.per.machine", "32")));	
	        for (int i = 0; i < regions.size(); i++)
	        {
	        	String[] region = regions.get(i);
	        	Runnable deleteWorker = new S3DeleteWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".final.bam", logPath);
				es.execute(deleteWorker);
		    	deleteWorker = new S3DeleteWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".final.bai", logPath);
				es.execute(deleteWorker);
		    	deleteWorker = new S3DeleteWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".final.bed", logPath);
		    	es.execute(deleteWorker);
	        	deleteWorker = new S3DeleteWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".raw.vcf.gz", logPath);
				es.execute(deleteWorker);
		    	deleteWorker = new S3DeleteWorker( loader, m_S3DownloadURL, m_fileName + "." + region[0] + ":" + region[1] + "-" + region[2] + ".raw.vcf.gz.tbi", logPath);
				es.execute(deleteWorker);
		    }
			
			es.shutdown();
			while (!es.isTerminated()) {}
	
	        System.out.println(new Timestamp(date.getTime()) + ": S3 delete is done.");
		}
	}
}