package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import java.util.Date;
import java.sql.Timestamp;

/**
 * GroupSampleCaller calls GroupVCFWorker to group individual VCF files into 
 * one family.
 * 
 * @author Guorong Xu
 */
public class GroupSampleCaller {
	/**
	 * Creates a shell command to call GroupVCFWorker to group individual VCF
	 * files into one family.
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
			
			HashMap<String, String> m_sampleGroups = parser.getSampleGroups();
			String S3UploadURL = parser.getS3UploadURL();
			
			if ( ! m_sampleGroups.isEmpty() )
			{
				for (String group : m_sampleGroups.keySet()) 
				{
					Date date= new Date();
					File dir = new File(rootPath + "results/" + group);
					if (!dir.exists())
						dir.mkdirs();
					
					String samples = m_sampleGroups.get(group);
					String bashFile = rootPath + "results/" + group + "/groupworker_" + group + ".sh";				
					FileWriter filewriter = new FileWriter(bashFile, false);
					
					//the same URL
					filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=groupworker S3UploadURL=" + S3UploadURL + "/" + group
							+ " S3DownloadURL=" + S3UploadURL  + " fileName1=" + group + " fileName2=" + samples + " description= date= configFile=" + homeDirctory + "/system.conf\n");
					
					/* Create a command for head node to send out jobs to work 
					 * nodes.
					 * default slot of group.call is 2*/
					String command = " qsub -pe smp " + loader.get("group.call.slot", "4") + " " + bashFile;
					System.out.println(new Timestamp(date.getTime()) + command);
					
					Runtime rt = Runtime.getRuntime();
					Process pr = rt.exec(command);
					
					filewriter.close();
				}
			}
			/*track the PBS queue statue*/
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "groupwork");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
