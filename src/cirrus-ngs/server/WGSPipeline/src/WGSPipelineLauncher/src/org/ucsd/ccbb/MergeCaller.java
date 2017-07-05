package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import java.util.Date;
import java.sql.Timestamp;
/**
 * MergeCaller calls MergeWorker to merge splitted final BAM  files and VCF 
 * files.
 * 
 * @author Guorong Xu
 */
public class MergeCaller {

	/**
	 * Create a shell command to call MergeWorker to merge splitted final BAM 
	 * files and VCF files.
	 * @param rootPath		a path points to the work space
	 * @param loader		load the configuration file
	 * @param yamlFile		a yaml File that includes demands of customer
	 */
	public void run(String rootPath, ConfigLoader loader, String yamlFile) 
	{
		String slots = loader.get("slots.per.machine", "32");
		/*get the WGS pipeline home directory*/
		String homeDirctory = loader.get("home.dirctory", "/shared/workspace/WGSPipeline");
		/*get the WGS pipeline manager*/
		String wgsmanager = loader.get("WGS.Pipeline.wgsmanager", "/shared/workspace/WGSPipeline/libs/WGSPipeline.jar");
			        
		try
		{		
			YamlParser parser = YamlParser.getInstance();
			
			String S3UploadURL = parser.getS3UploadURL();
			ArrayList<String[]> fastqList = parser.getFastqList();
			
		    for ( int index = 0; index < fastqList.size(); index++ ) 
		    {
				Date date= new Date();
				String fileName = fastqList.get(index)[0];
				String S3DownloadURL = fastqList.get(index)[1];
				File dir = new File(rootPath + "results/" + fileName);
				if (!dir.exists())
					dir.mkdirs();
				
				String bashFile = rootPath + "results/" + fileName + "/mergeworker_" + fileName + ".sh";
				FileWriter filewriter = new FileWriter(bashFile, false);
				
				/*Create a command*/
				filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=mergeworker S3UploadURL=" + S3UploadURL + "/" + fileName 
						+ " S3DownloadURL=" + S3UploadURL + "/" + fileName + " fileName1=" + fileName + " fileName2= description= date= configFile=" 
						+ homeDirctory + "/system.conf\n");	
				
				/* Create a command for head node to send out jobs to work 
				 * nodes.
				 */
				String command = " qsub -pe smp " + slots + " " + bashFile;
				System.out.println(new Timestamp(date.getTime()) + command);
				
				Runtime rt = Runtime.getRuntime();
				/*Executes the command specified before in the process pr*/
				Process pr = rt.exec(command);
				
				filewriter.close();	
		    }
		    
		    /*track the PBS queue statue*/
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "mergework");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
