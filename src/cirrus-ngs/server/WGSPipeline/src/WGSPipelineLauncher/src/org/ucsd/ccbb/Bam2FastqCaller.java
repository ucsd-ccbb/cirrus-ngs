package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;

/**
 * Bam2FastqCaller creates a command to convert a raw BAM file to paired-end 
 * fastq files.
 * 
 * @author Guorong Xu
 */
public class Bam2FastqCaller {

	/**
	 * Converts a raw BAM file to a paired-end fastq file
	 * @param rootPath  a path points to the work space
	 * @param loader	load the configuration file
	 * @param yamlFile	a yaml File that includes demands of customer
	 * @throws Exception
	 */
	public void run(String rootPath, ConfigLoader loader, String yamlFile) throws Exception
	{
		/*get the WGS pipeline home directory*/
		String homeDirctory = loader.get("home.dirctory", "/scratch/workspace");
		/*get the WGS pipeline manager*/
		String wgsmanager = loader.get("WGS.Pipeline.wgsmanager", "/shared/workspace/WGSPipeline/libs/WGSPipeline.jar");	
		
		try
		{	
			/* get the project date, S3 upload URL and BAM files informations 
			 * from yaml file
			 */
			YamlParser parser = YamlParser.getInstance();
			String projectDate = parser.getProjectDate();
			String S3UploadURL = parser.getS3UploadURL();
			ArrayList<String[]> bamList = parser.getBamList();
			
		    for ( int index = 0; index < bamList.size(); index++ ) 
		    {
		    		Date date= new Date();
				String fileName = bamList.get(index)[0];
				String S3DownloadURL = bamList.get(index)[1];
				String description = bamList.get(index)[2];
				File dir = new File(rootPath + "results/" + fileName);
				if (!dir.exists())
					dir.mkdirs();
				
				/* Create a bash file
				 * Data will be written to the beginning of the file
				 */
				String bashFile = rootPath + "results/" + fileName + "/bam2fastqworker_" + fileName + ".sh";
				
				FileWriter filewriter = new FileWriter(bashFile, false);
			
				String bamFile = fileName + ".bam";
				
				/*Create a command to covert bam file to fastq file in the bash file*/
				filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=bam2fastq S3UploadURL=" + S3UploadURL
						+ " S3DownloadURL=" + S3DownloadURL + " fileName1=" + bamFile + " fileName2= description=" + description 
						+ " date=" + projectDate + " configFile=" + homeDirctory + "/system.conf\n");	
				
				/* Create a command for head node to send out jobs to work nodes.
				 * Default slot for bam2fastq is 4
				 */ 
				String command = " qsub -pe smp " + loader.get("bam2fastq.slot", "4") + " " + bashFile;
				System.out.println(new Timestamp(date.getTime()) + command);
				
				/*Get the current runtime*/
				Runtime rt = Runtime.getRuntime();
				/*Executes the command in the process pr*/
				Process pr = rt.exec(command);
				
				filewriter.close();	
		    }
		    
		    /*track the PBS queue statue*/
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "bam2fastq");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
