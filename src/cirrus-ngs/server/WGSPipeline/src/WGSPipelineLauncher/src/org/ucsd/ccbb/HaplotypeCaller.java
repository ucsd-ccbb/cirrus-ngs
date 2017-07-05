package org.ucsd.ccbb;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import java.util.Date;
import java.sql.Timestamp;
/**
 * HaplotypeCaller calls GATKHaplotypeWorker to apply HaplotypeCaller in GATK
 * 
 * @author Guorong Xu
 */
public class HaplotypeCaller {
	
	/**an array of string arrays stores 25 chromosome with their names and 
	 * corresponding lengthes.*/
	String[][] m_chromLength = new String[][]{		
			{"chr1", "249250621"},
			{"chr2", "243199373"},
			{"chr3", "198022430"},
			{"chr4", "191154276"},
			{"chr5", "180915260"},
			{"chr6", "171115067"},
			{"chr7", "159138663"},
			{"chr8", "146364022"},
			{"chr9", "141213431"},
			{"chr10", "135534747"},
			{"chr11", "135006516"},
			{"chr12", "133851895"},
			{"chr13", "115169878"},
			{"chr14", "107349540"},
			{"chr15", "102531392"},
			{"chr16", "90354753"},
			{"chr17", "81195210"},
			{"chr18", "78077248"},
			{"chr19", "59128983"},
			{"chr20", "63025520"},
			{"chr21", "48129895"},
			{"chr22", "51304566"},
			{"chrX", "155270560"},
			{"chrY", "59373566"},
			{"chrM", "16571"},
	};
	
	/**
	 * Creates a shell command to call GATKHaplotypeWorker to apply 
	 * HaplotypeCaller in GATK.
	 * @param rootPath	a path points to the work space
	 * @param loader	load the configuration file
	 * @param yamlFile	a yaml File that includes demands of customer
	 */
	public void run(String rootPath, ConfigLoader loader, String yamlFile) 
	{	
		/*get the WGS pipeline home directory*/
		String homeDirctory = loader.get("home.dirctory", "/shared/workspace/WGSPipeline");
		/*get the WGS pipeline manager*/ 
		String wgsmanager = loader.get("WGS.Pipeline.wgsmanager", "/shared/workspace/WGSPipeline/libs/WGSPipeline.jar");
		/*get the sam.chrom.length.interval, or assign 30000000 to the field interval*/
		int interval = Integer.valueOf(loader.get("sam.chrom.length.interval", "30000000"));
		
		try
		{		
			/* get the S3 upload URL and fastq files informations from yaml 
			 * file
			 */
			YamlParser parser = YamlParser.getInstance();
			
			String S3UploadURL = parser.getS3UploadURL();
			ArrayList<String[]> fastqList = parser.getFastqList();
			
		    for ( int index = 0; index < fastqList.size(); index++ ) 
		    {
				String fileName = fastqList.get(index)[0];
				String group = fastqList.get(index)[3];
				
				File dir = new File(rootPath + "results/" + fileName);
				if (!dir.exists())
					dir.mkdirs();
				
				//m_chromLength.length is 25
				for (int i = 0; i < m_chromLength.length; i ++)
				{
					String chrom = m_chromLength[i][0];
					/* get the integer representation of the second string in 
					 * each element of the string array
					 */
					int length = Integer.valueOf(m_chromLength[i][1]);
					int regionNum = length / interval + 1;
					
					for (int j = 0; j < regionNum; j ++)
					{
						Date date= new Date();
						//start with 1
						int start = j * interval + 1;
						int end = (j + 1) * interval;
						
						if (end > length)
							end = length;
							
						/*Repetitive variable name*/
						String splitBAMFile = fileName + "."  + chrom + ":" + start + "-" + end;	
						String bashFileName = fileName + "."  + chrom + "." + start + "-" + end;
						/* Create a bash file
						 * Data will be written to the beginning of the file
						 */
						String bashFile = rootPath + "results/" + fileName + "/haplotypeworker_" + bashFileName + ".sh";
						FileWriter filewriter = new FileWriter(bashFile, false);
						
						//Once we uploaded the split BAM files to S3, then we should always download the Intermediate files using S3UploadURL.
						filewriter.write("#!/bin/bash\n");	
						filewriter.write("java -Xms454m -Xmx2g -jar " + wgsmanager + " command=haplotypeworker S3UploadURL=" + S3UploadURL + "/" + fileName 
								+ " S3DownloadURL=" + S3UploadURL + "/" + fileName + " fileName1=" + splitBAMFile + ".final" 
								+ " fileName2=" + group + " description= date= configFile=" + homeDirctory + "/system.conf\n");	
								
						filewriter.close();
						
						/* Create a command for head node to send out jobs to work 
						 * nodes.
						 * default slot of haplotype.call is 2*/
						String command = " qsub -pe smp " + loader.get("haplotype.call.slot", "2") + " " + bashFile;
						if (!"NA".equalsIgnoreCase(group))
							command = " qsub -pe smp " + loader.get("haplotype.group.call.slot", "4") + " " + bashFile;
						
						System.out.println(new Timestamp(date.getTime()) + command);
						
						/*Get the current runtime*/
						Runtime rt = Runtime.getRuntime();
						/*Executes the command specified before in the process pr*/
						Process pr = rt.exec(command);
					}
				}				
		    }
		    
		    /*track the PBS queue statue*/
		    new PBSTracker().trackPBSQueue(Integer.valueOf(loader.get("trackTime", "1")), "haplotype");	
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
