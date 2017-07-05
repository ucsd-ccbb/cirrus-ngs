package org.ucsd.ccbb;

import java.sql.Timestamp;
import java.util.Date;

/**
 * WGSManager processes command to call corresponding workers to finish 
 * different commands. The commands include 
 * <ul>
 * <li>converting BAM to Fastq
 * <li>calling fastqc?
 * <li>BWAAlignmentWorker
 * <li>PostAlignmentWorker
 * <li>GATKHaplotypeWorker
 * <li>MergeWorker 
 * <li>GroupVCFWorker
 * <li>VariantFilterWorker.
 * </ul>
 * @author Guorong Xu
 */
public class WGSManager {
	/**a string stores the command type.
	 * The command type includes 
	 * <ul>
	 * <li>"bam2fastq"
	 * <li>"fastqcworker"
	 * <li>"alignworker"
	 * <li>"postalignworker"
	 * <li>"haplotypeworker"
	 * <li>"mergeworker"
	 * <li>"groupworker"  
	 * <li>"filterworker"
	 * </ul>
	 */ 
	public String m_commandType = "";
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores the first file name of the paired-end fastq file or a 
	 * BAM file name
	 */
	public String m_fileName1 = "";
	/**a string stores the second file name of the paired-end fastq file or the 
	 * description
	 */
	public String m_fileName2 = "";
	/**a string stores the configuration file?*/
	public String m_configFile = "";
	/**a string stores the description of a sample*/
	public String m_description = "";
	/**a string stores the date of the project*/
	public String m_date = "";
	
	/**
	 * The main method constructs a WGSManager instance based on 8 input 
	 * parameters.
	 * @param args  The length of the string array is 8. 
	 * 				It includes the command type, the S3 upload URL, the S3 
	 * 				download URL, fileName1, fileName2, the description of
	 * 				a sample, the date of the project and the configuration
	 * 				file 
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception 
	{	
		for (int i = 0; i < args.length; i++)
		{
			System.out.println("args[" + i + "]: " + args[i]);
		}
		String commandType = args[0].substring(args[0].indexOf("=") + 1);
		String S3UploadURL = args[1].substring(args[1].indexOf("=") + 1);
		String S3DownloadURL = args[2].substring(args[2].indexOf("=") + 1);
		String fileName1 = args[3].substring(args[3].indexOf("=") + 1);
		String fileName2 = args[4].substring(args[4].indexOf("=") + 1);
		String description = args[5].substring(args[5].indexOf("=") + 1);
		String date = args[6].substring(args[6].indexOf("=") + 1);
		String configFile = args[7].substring(args[7].indexOf("=") + 1);
			
		WGSManager worker = new WGSManager(commandType, S3UploadURL, S3DownloadURL, 
				fileName1, fileName2, description, date, configFile);
		worker.processCommand();
	}
	
	/**
	 * Construct a new WGSManager with specified commandType, S3UploadURL, 
	 * S3DownloadURL, fileName1, fileName2, description, date and configFile
	 * @param commandType	The command type
	 * @param S3UploadURL	The S3 upload URL
	 * @param S3DownloadURL	The S3 download URL
	 * @param fileName1 	the first file name of the paired-end fastq file or
	 * 						a BAM file name
	 * @param fileName2		the second file name for the paired- end fastq file 
	 * 						or the description
	 * @param description	The description of a sample
	 * @param date			The date of the project
	 * @param configFile	The path of the configuration file
	 */
	public WGSManager(String commandType, String S3UploadURL, String S3DownloadURL, 
			String fileName1, String fileName2, String description, String date, String configFile) {
		this.m_commandType = commandType;
		this.m_S3UploadURL = S3UploadURL;
		this.m_S3DownloadURL = S3DownloadURL;
		this.m_fileName1 = fileName1;
		this.m_fileName2 = fileName2;
		this.m_description = description;
		this.m_date = date;
		this.m_configFile = configFile;
	}

	/**
	 * According to the command type, this method creates corresponding worker 
	 * instance and process command further. 
	 *  
	 * @throws Exception
	 */
	public void processCommand() throws Exception {
		Date date= new Date();
		System.out.println(new Timestamp(date.getTime()) + "Command type: " + m_commandType + ".");
		switch (m_commandType) {
		/*If the command type is bam2fastq, create a BAM2FastqWorker instance*/
		case "bam2fastq":
			BAM2FastqWorker b2fworker = new BAM2FastqWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_configFile);
			b2fworker.processCommand();
			break;
		/*If the command type is FastqcWorker, create a FastqcWorker instance*/	
		case "fastqcworker":
			FastqcWorker qcworker = new FastqcWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_fileName2, 
					m_description, m_date, m_configFile);
			qcworker.processCommand();
			break;
		/*If the command type is alignworker, create a BAMAlignmentWorker instance*/
			//From this point, it is a separate part?*/
		case "alignworker":
			BWAAlignmentWorker fworker = new BWAAlignmentWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_fileName2, 
					m_description, m_date, m_configFile);
			fworker.processCommand();
			break;
		/*If the command type is postalignworker, create a PostAlignmentWorker instance*/
		case "postalignworker":
			PostAlignmentWorker bworker = new PostAlignmentWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_fileName2, m_configFile);
			bworker.processCommand();
			break;
		/*If the command type is haplotypeworker, create a GATKHaplotypeWorker instance*/
		case "haplotypeworker":
			GATKHaplotypeWorker gatkworker = new GATKHaplotypeWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_fileName2, m_configFile);
			gatkworker.processCommand();
			break;
		/*If the command type is mergeworker, create a MergeWorker instance*/
		case "mergeworker":			
			MergeWorker mworker = new MergeWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_configFile);
			mworker.processCommand();
			break;
		/*If the command type is groupworker, create a GroupVCFWorker instance*/
		case "groupworker":			
			GroupVCFWorker gworker = new GroupVCFWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_fileName2, 
					m_configFile);
			gworker.processCommand();
			break;
		/*If the command type is filterworker, create a VariantFilterWorker instance*/
		case "filterworker":			
			VariantFilterWorker vfworker = new VariantFilterWorker(m_S3UploadURL, m_S3DownloadURL, m_fileName1, m_configFile);
			vfworker.processCommand();
			break;
		/*If the command type is none of above, print out "INVALID COMMAND"*/
		default:
			System.out.println("INVALID COMMAND!");
		}
	}
}