package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * GroupVCFWorker groups individual VCF files into one family.
 * 
 * @author Guorong Xu
 */
public class GroupVCFWorker {
	
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	public String m_groupFileName = "";
	/**a string stores a list of VCF files*/
	public String m_vcfFileList = "";
	/**a string stores the configuration file name?*/
	public String m_configFile = "";
	
	/**
	 * Constructs a GroupVCFWorker instance with specified S3 upload URL, S3 
	 * download URL, group name, a list of VCF files and configuration file.
	 * @param S3UploadURL		the s3 upload URL
	 * @param S3DownloadURL		the s3 download URL
	 * @param groupName			the group name?
	 * @param vcfFileList		a list of VCF files
	 * @param configFile		the configuration file
	 */
	public GroupVCFWorker(String S3UploadURL, String S3DownloadURL, String groupName, String vcfFileList, String configFile) 
	{
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_groupFileName = groupName;
		this.m_vcfFileList = vcfFileList;
		this.m_configFile = configFile;
	}
	
	/**
	 * Processes commands.
	 * First, download files from S3. Then, group samples into one family. Next,
	 * upload files to S3. Last, delete VCF files from workspace.
	 * @throws Exception
	 */
	public void processCommand() throws Exception {		
		Date date= new Date();
		
		System.out.println(new Timestamp(date.getTime()) + ": System is processing " + m_groupFileName + "," + m_vcfFileList);
		ConfigLoader loader = new ConfigLoader();
		loader.parsePathFile(m_configFile);
		
		String workspace = loader.get("work.space", "/scratch/workspace");	
		String homeDirctory = loader.get("home.dirctory", "/shared/workspace/WGSPipeline");
		String logPath = homeDirctory + "/results/logs";		
		File dir = new File(logPath);
		if (!dir.exists())
			dir.mkdirs();
		
		ArrayList<String> vcfFiles = getVCFList();
        ArrayList<String> vcfFileList = new ArrayList<String>();
        
        /*download files from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
	    ExecutorService es = Executors.newFixedThreadPool(vcfFiles.size());	
        for (int i = 0; i < vcfFiles.size(); i++)
        {
        	vcfFileList.add(vcfFiles.get(i) + ".final.vcf.gz ");
 			Runnable worker = new S3DownloadWorker( loader, m_S3DownloadURL + "/" + vcfFiles.get(i), vcfFiles.get(i) + ".final.vcf.gz", logPath);
			es.execute(worker);
	    }		
		es.shutdown();
		while (!es.isTerminated()) {}
		
		es = Executors.newFixedThreadPool(vcfFiles.size());	
        for (int i = 0; i < vcfFiles.size(); i++)
        {
 			Runnable worker = new S3DownloadWorker( loader, m_S3DownloadURL + "/" + vcfFiles.get(i), vcfFiles.get(i) + ".final.vcf.gz.tbi", logPath);
			es.execute(worker);
	    }		
		es.shutdown();
		while (!es.isTerminated()) {}
        
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
        
		// Grouping samples into one family
		System.out.println(new Timestamp(date.getTime()) + ": Grouping VCF is processing.");
		GroupVCFExecutor gvcf = new GroupVCFExecutor();
		gvcf.setShellFile(loader.get("group.vcf.shell", "/shared/workspace/WGSPipeline/scripts/group_vcf.sh"));			
		gvcf.setGroupName(m_groupFileName);
		gvcf.setVcfList(vcfFileList);
		gvcf.setLogFile(logPath);
		gvcf.execute("", "");
        System.out.println(new Timestamp(date.getTime()) + ": Grouping VCF is done.");
		        
        /*upload files to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");
	    S3UploadExecutor s3u= new S3UploadExecutor();
	    s3u.setShellFile(loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setS3UploadURL(m_S3UploadURL);
	    s3u.setLogFile(logPath);	    
	    s3u.setFileName(m_groupFileName + ".g.vcf.gz");
	    s3u.execute("", m_groupFileName + ".g.vcf.gz");
	    
	    s3u.setFileName(m_groupFileName + ".g.vcf.gz.tbi");
	    s3u.execute("", m_groupFileName + ".g.vcf.gz.tbi");
        System.out.println(new Timestamp(date.getTime()) + ": S3 upload is done.");
        
		//Deleting VCF files from workspace
	    System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is processing.");
	    DeleteWorkspaceExecutor deleter = new DeleteWorkspaceExecutor();
	    deleter.setShellFile(loader.get("delete.workspace.shell", "/shared/workspace/WGSPipeline/scripts/delete.sh"));
        for (int i = 0; i < vcfFiles.size(); i++)
        {    
		    deleter.setInputFile(workspace + "/" + vcfFiles.get(i));
		    deleter.execute("", workspace + "/" + vcfFiles.get(i));
        }
		System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is done.");   
	}
	
	/**
	 * Separate a list of VCF files into individual VCF files.
	 * @return	an arraylist of VCF files. Each element of the arraylist is a 
	 * 			VCF file
	 */
	public ArrayList<String> getVCFList()
	{
		ArrayList<String> vcfList = new ArrayList<String>();
		
		int start = 0;
		int next = m_vcfFileList.indexOf("/");
		while(next >= 0) {
			// "fileName1/fileName2/"
		     String fileName = m_vcfFileList.substring(start, next);
		     vcfList.add(fileName);
		     
		     start = next + 1;
		     next = m_vcfFileList.indexOf("/", start);
		     
		}
		
		return vcfList;
	}
}