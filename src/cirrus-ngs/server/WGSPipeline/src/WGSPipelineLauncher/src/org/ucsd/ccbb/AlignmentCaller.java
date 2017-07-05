package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Date;
import java.sql.Timestamp;

/**
 * AlignmentCaller calls BWA to map fastq to human reference genome(hg19)
 * 
 * @author Guorong Xu
 */
public class AlignmentCaller {
	public HashMap<String, ArrayList<String>> m_sampleGroups = new HashMap<String, ArrayList<String>>();
	
	/**
	 * Creates a shell command to call BWAAlignmentWorker to process all 
	 * alignments, sorting and splitting
	 * @param rootPath  a path points to the workspace
	 * @param loader	load the configuration file
	 * @param yamlFile	a yaml File that includes demands of customer
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
			/* get the project date, S3 upload URL and fastq files informations 
			 * from yaml file
			 */
			YamlParser parser = YamlParser.getInstance();
			String projectDate = parser.getProjectDate();
			String S3UploadURL = parser.getS3UploadURL();
			ArrayList<String[]> fastqList = parser.getFastqList();
			
		    for ( int index = 0; index < fastqList.size(); index++ ) 
		    {
		    		Date date= new Date();
				String fileName = fastqList.get(index)[0];
				String S3DownloadURL = fastqList.get(index)[1];
				String description = fastqList.get(index)[2];
				File dir = new File(rootPath + "results/" + fileName);
				if (!dir.exists())
					dir.mkdirs();
				
				/* Create a bash file
				 * Data will be written to the beginning of the file
				 */
				String bashFile = rootPath + "results/" + fileName + "/alignworker_" + fileName + ".sh";
				FileWriter filewriter = new FileWriter(bashFile, false);
			
				String fastq1 = fileName + "_R1.fq.gz";
				String fastq2 = fileName + "_R2.fq.gz";
				
				filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=alignworker S3UploadURL=" + S3UploadURL
						+ " S3DownloadURL=" + S3DownloadURL + " fileName1=" + fastq1 + " fileName2=" + fastq2 
						+ " description=" + description + " date=" + projectDate + " configFile=" + homeDirctory + "/system.conf\n");	
				
				String command = " qsub -pe smp " + slots + " " + bashFile;
				System.out.println(new Timestamp(date.getTime()) + command);
				
				Thread.sleep(5000);
				
				/*Get the current runtime*/
				Runtime rt = Runtime.getRuntime();
				/*Executes the command specified before in the process pr*/
				Process pr = rt.exec(command);
				
				filewriter.close();	
		    }
		    /*track the PBS queue statue*/
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "alignwork");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
