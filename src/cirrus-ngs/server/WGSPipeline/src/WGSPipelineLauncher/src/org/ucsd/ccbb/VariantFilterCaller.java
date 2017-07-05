package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import java.util.Date;
import java.sql.Timestamp;
/**
 * VariantFilterCaller uses haplotype to call variants
 * 
 * @author Guorong Xu
 */
public class VariantFilterCaller {
	
	/**
	 * Creates a shell command to call VariantFilterWorker to filter variant
	 * using GATK VQSR filteration.
	 * @param rootPath		a path points to the work space
	 * @param loader		load the configuration file
	 * @param yamlFile		a yaml File that includes demands of customer 
	 */
	public void run(String rootPath, ConfigLoader loader, String yamlFile)
	{
		try
		{		
			/*get the WGS pipeline home directory*/
			String homeDirctory = loader.get("home.dirctory", "/shared/workspace/WGSPipeline");	
			/*get the WGS pipeline manager*/
			String wgsmanager = loader.get("WGS.Pipeline.wgsmanager", "/shared/workspace/WGSPipeline/libs/WGSPipeline.jar");	
			
			YamlParser parser = YamlParser.getInstance();
			
			String S3UploadURL = parser.getS3UploadURL();
			HashMap<String, String> sampleGroups = parser.getSampleGroups();
			ArrayList<String[]> fastqList = parser.getFastqList();
			
			/*if the sample group is not empty*/
			if ( ! sampleGroups.isEmpty() )
			{
				for (String group : sampleGroups.keySet()) 
				{
					Date date= new Date();
					File dir = new File(rootPath + "results/" + group);
					if (!dir.exists())
						dir.mkdirs();
					String bashFile = rootPath + "results/" + group + "/filterworker_" + group + ".sh";					
					FileWriter filewriter = new FileWriter(bashFile, false);
					
					filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=filterworker S3UploadURL=" + S3UploadURL + "/" + group
							+ " S3DownloadURL=" + S3UploadURL + "/" + group +" fileName1=" + group + ".g" + " fileName2= description= date= configFile=" 
							+ homeDirctory + "/system.conf\n");	
					
					/* Create a command for head node to send out jobs to work nodes.
					 * Default slot of variant.filter is 4
					 */ 
					String command = " qsub -pe smp " + loader.get("variant.filter.slot", "4") + " " + bashFile;
					System.out.println(new Timestamp(date.getTime()) + command);
					
					Runtime rt = Runtime.getRuntime();
					/*Executes the command in the process pr*/
					Process pr = rt.exec(command);
					
					filewriter.close();
				}
			}
			/*if the sample group is empty*/
			else
			{
			    for ( int index = 0; index < fastqList.size(); index++ ) 
			    {
					Date date2= new Date();
					String fileName = fastqList.get(index)[0];	
					File dir = new File(rootPath + "results/" + fileName);
					if (!dir.exists())
						dir.mkdirs();
					
					String bashFile = rootPath + "results/" + fileName + "/filterworker_" + fileName + ".sh";
					FileWriter filewriter = new FileWriter(bashFile, false);
									
					filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=filterworker S3UploadURL=" + S3UploadURL + "/" + fileName 
							+ " S3DownloadURL=" + S3UploadURL + "/" + fileName + " fileName1=" + fileName + ".final" + " fileName2= description= date= configFile=" 
							+ homeDirctory + "/system.conf\n");
					
					/* Create a command for head node to send out jobs to work nodes.
					 * Default slot of variant.filter is 4
					 */ 
					String command = " qsub -pe smp " + loader.get("variant.filter.slot", "4") + " " + bashFile;
					System.out.println(new Timestamp(date2.getTime()) + command);
					
					Runtime rt = Runtime.getRuntime();
					Process pr = rt.exec(command);
					
					filewriter.close();
			    }
			}
			
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "filterwo");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
