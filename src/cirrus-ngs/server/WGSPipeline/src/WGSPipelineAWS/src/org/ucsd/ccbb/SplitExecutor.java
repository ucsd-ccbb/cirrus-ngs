package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * SplitExecutor splits a BWA file to small BAM files in the specific region.
 * 
 * @author Guorong Xu
 */
public class SplitExecutor extends CommandExecutor {
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the SAM file?*/
	public String m_samFile = "";
	/**a string stores the name of a chromosome*/
	public String m_chrom = "";
	/**an integer stores the start point of a subchromosome*/
	public int m_start = 0;
	/**an integer stores the end point of a subchromosome*/
	public int m_end = 0;
	/**a string stores the log path*/
	public String m_logFile = "";
	
	/**
	 * This function override the function assembleCommands() from class 
	 * CommandExecutor. 
	 * This function assembles a shell command with shell script path?
	 * @param filePath	unused in this function
	 * @param fileName	unused in this function
	 * @return	a string array includes the command
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(m_samFile);
		commandArray.add(m_chrom);
		commandArray.add(String.valueOf(m_start));
		commandArray.add(String.valueOf(m_end));
		commandArray.add(m_logFile);

		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * This function set the input path as the shell script path.
	 * @param fileName	The path will be set to the field m-shellScript
	 */
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	/**
	 * Stores the input SAM file name.
	 * @param samFile	The SAM file name will be set to the field m_samFile
	 */
	public void setSAMFile(String samFile)
	{
		m_samFile = samFile;
	}
	
	/**
	 * This function specifies a region with the name of a chromosome, the 
	 * start point of a subchromosome on the chromosome and the end point of 
	 * the subchromosome.
	 * @param chrom		the name of a chromosome
	 * @param start		the start point of a subchromosome on the chromosome
	 * @param end		the end point of a subchromosome on the chromosome
	 */
	public void setRegion(String chrom, int start, int end)
	{
		m_chrom = chrom;
		m_start = start;
		m_end	= end;
	}
	
	/**
	 * This function set the input path as the log path.
	 * @param logFile	The path will be set to the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
}
