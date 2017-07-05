package org.ucsd.ccbb;

import java.io.File;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * VariantFilterWorker filters variant using GATK VQSR filteration.
 * 
 * @author Guorong Xu
 */
public class VariantFilterWorker {
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores a list of VCF files?*/
	public String m_vcfFileName = "";
	/**a string stores the configuration file name?*/
	public String m_configFile = "";
	
	/**
	 * Constructs a VariantFilterWorker instance with specified S3 upload URL, S3 
	 * download URL, VCF file name and configuration file.
	 * @param S3UploadURL
	 * @param S3DownloadURL
	 * @param vcfFileName
	 * @param configFile
	 */
	public VariantFilterWorker(String S3UploadURL, String S3DownloadURL, String vcfFileName, String configFile) 
	{
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_vcfFileName = vcfFileName;
		this.m_configFile = configFile;
	}
	
	/**
	 * Processes commands
	 * First, download files from S3. Then, filter variants. Next, upload files
	 * to S3. Last, delete VCF files from workspace.
	 * @throws Exception
	 */
	public void processCommand() throws Exception {		
		Date date= new Date();
		
		System.out.println(new Timestamp(date.getTime()) + ": System is processing " + m_vcfFileName + ".");
		ConfigLoader loader = new ConfigLoader();
		loader.parsePathFile(m_configFile);
		
		String workspace = loader.get("work.space", "/scratch/workspace");	
		String homeDirctory = loader.get("home.dirctory", "/shared/workspace/WGSPipeline");
		String logPath = homeDirctory + "/results/logs";		
		File dir = new File(logPath);
		if (!dir.exists())
			dir.mkdirs();
        
		/*download files from S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 download is processing.");
	    S3DownloadExecutor s3d= new S3DownloadExecutor();
	    s3d.setShellFile(loader.get("s3download.shell", "/shared/workspace/WGSPipeline/scripts/s3download.sh"));
	    s3d.setS3DownloadURL(m_S3DownloadURL);
	    s3d.setLogFile(logPath);
	    s3d.setFileName(m_vcfFileName + ".vcf.gz");
	    s3d.execute("", m_vcfFileName + ".vcf.gz");
        
	    s3d.setFileName(m_vcfFileName + ".vcf.gz.tbi");
	    s3d.execute("", m_vcfFileName + ".vcf.gz.tbi");
        System.out.println(new Timestamp(date.getTime()) + ": S3 download is done.");
        
		/*filter variants*/
		System.out.println(new Timestamp(date.getTime()) + ": VariantFilter is processing.");
		VariantFilterExecutor vfe = new VariantFilterExecutor();
		vfe.setShellFile(loader.get("filter.variant.shell", "/shared/workspace/WGSPipeline/scripts/filter_variant.sh"));			
		vfe.setFileName(m_vcfFileName);
		vfe.setLogFile(logPath);
		vfe.execute("", m_vcfFileName);
        System.out.println(new Timestamp(date.getTime()) + ": VariantFilter is done.");
		      
        /*upload files to S3*/
	    System.out.println(new Timestamp(date.getTime()) + ": S3 upload is processing.");
	    S3UploadExecutor s3u= new S3UploadExecutor();
	    s3u.setShellFile(loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setS3UploadURL(m_S3UploadURL);
	    s3u.setLogFile(logPath);	    
	    s3u.setFileName(m_vcfFileName + ".vqsr.vcf");
	    s3u.execute("", m_vcfFileName + ".vqsr.vcf");
	    
	    s3u.setFileName(m_vcfFileName + ".vqsr.vcf.tbi");
	    s3u.execute("", m_vcfFileName + ".vqsr.vcf.tbi");
        System.out.println(new Timestamp(date.getTime()) + ": S3 upload is done.");
        
		//Deleting VCF files from workspace
	    System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is processing.");
	    DeleteWorkspaceExecutor deleter = new DeleteWorkspaceExecutor();
	    deleter.setShellFile(loader.get("delete.workspace.shell", "/shared/workspace/WGSPipeline/scripts/delete.sh")); 
	    deleter.setInputFile(workspace + "/" + m_vcfFileName);
	    deleter.execute("", workspace + "/" + m_vcfFileName);
		System.out.println(new Timestamp(date.getTime()) + ": Delete-Workspace is done.");      

	}
}