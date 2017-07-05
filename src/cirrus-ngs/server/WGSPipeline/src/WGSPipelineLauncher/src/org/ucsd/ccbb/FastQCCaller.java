package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import java.util.Arrays;
import java.util.Date;
import java.sql.Timestamp;
/**
 * FastQCCaller calls FastqWorker to run fastQC
 * 
 * @author Guorong Xu
 */
public class FastQCCaller {
	
	/**
	 * Creates a shell command to call FastqWorker to run fastQC
	 * @param rootPath 	a path points to the work space
	 * @param loader	load the configuration file
	 * @param yamlFile  a yaml File that includes demands of customer
	 */
	public void run(String rootPath, ConfigLoader loader, String yamlFile) 
	{	
		/*get the WGS pipeline home directory*/
		String homeDirctory = loader.get("home.dirctory", "/shared/workspace/WGSPipeline");
		/*get the WGS pipeline manager directory*/
		String wgsmanager = loader.get("WGS.Pipeline.wgsmanager", "/shared/workspace/WGSPipeline/libs/WGSPipeline.jar");		
		
		try
		{	
			/* get the project date, S3 upload URL and fastq files informations 
			 * from yaml file
			 */
			Date date= new Date();
			YamlParser parser = YamlParser.getInstance();
			String projectDate = parser.getProjectDate();
			String S3UploadURL = parser.getS3UploadURL();
			ArrayList<String[]> fastqList = parser.getFastqList();
			
			//---------------------------------------
			System.out.println(new Timestamp(date.getTime()) + " fastqList: ");
            for (int i=0; i<fastqList.size();i++){
                String[] temp=fastqList.get(i);
                System.out.println(Arrays.toString(temp));
            }
            //--------------------------------------------
			//System.out.println("fastqList: " + fastqList);
			
		    for ( int index = 0; index < fastqList.size(); index++ ) 
		    {
		    	Date date2= new Date();
				String fileName = fastqList.get(index)[0];
	    			
				System.out.println(fileName);
	    			
				String S3DownloadURL = fastqList.get(index)[1];
				String description = fastqList.get(index)[2];
				File dir = new File(rootPath + "results/" + fileName);
				if (!dir.exists())
					dir.mkdirs();
				
				/* Create a bash file
				 * Data will be written to the beginning of the file
				 */
				String bashFile = rootPath + "results/" + fileName + "/fastqc_" + fileName + ".sh";
				FileWriter filewriter = new FileWriter(bashFile, false);
			
				String fastq1 = fileName + "_R1.fq.gz";
				String fastq2 = fileName + "_R2.fq.gz";
				
				/*write a command into the bash file*/
				filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=fastqcworker S3UploadURL=" + S3UploadURL
						+ " S3DownloadURL=" + S3DownloadURL + " fileName1=" + fastq1 + " fileName2=" + fastq2 
						+ " description=" + description + " date=" + projectDate + " configFile=" + homeDirctory + "/system.conf\n");	
				
				/* Create a command for head node to send out jobs to work 
				 * nodes.
				 * default slot of fastqx.call is 2.
				 */
				String command = " qsub -pe smp " + loader.get("fastqc.call.slot", "2") + " " + bashFile;
				System.out.println(new Timestamp(date2.getTime()) + command);
				
				/*Get the current runtime*/
				Runtime rt = Runtime.getRuntime();
				/*Executes the command specified before in the process pr*/
				Process pr = rt.exec(command);
				
				filewriter.close();	
		    }
		    
		    /*track the PBS queue statue*/
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "fastqc");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
