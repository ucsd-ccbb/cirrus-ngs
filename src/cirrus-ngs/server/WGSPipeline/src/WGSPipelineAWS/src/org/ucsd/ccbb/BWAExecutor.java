package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * BWAExecutor executes BWA to run alignment for pair-end fastq files.
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class BWAExecutor extends CommandExecutor {
	
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the date of the project*/
	public String m_date = "";
	/**a string stores the description of the project*/
	public String m_description = "";
	/**a string stores the name of the first fastqFile path*/
	public String m_fastqFile1 = "";
	/**a string stores the name of the second fastqFile path*/
	public String m_fastqFile2 = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a string stores the suffix that will be searched in a file name*/
	public String m_suffix = "";
	
	/**
	 * The function assembles a shell command with shell script path, date, 
	 * description, the first fastq file path, the second fastq file path,
	 * file name and log file path
	 * @param filePath	unused in the function
	 * @param fileName	a file name that will be added to the command
	 * @return a string array includes the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		
		commandArray.add(m_shellScript);
		commandArray.add(m_date);
		commandArray.add(m_description);
		commandArray.add(m_fastqFile1);
		commandArray.add(m_fastqFile2);
		commandArray.add(fileName);
		commandArray.add(m_logFile);

		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * Stores the input path as the shell script path.
	 * @param fileName	The path will be set to the field m-shellScript
	 */
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	/**
	 * Stores the input string as the suffix 
	 * @param suffix	The string will be set to the field m_suffix
	 */
	public void setSuffix(String suffix)
	{
		m_suffix = suffix;
	}
	
	/**
	 * Stores the input string as the date of the project
	 * @param date		The string will be set to the field m-date
	 */
	public void setDate(String date)
	{
		m_date = date;
	}
	
	/**
	 * Stores the input string as the description
	 * @param description	The string will be set to the field m-description
	 */
	public void setDescription(String description)
	{
		m_description = description;
	}
	
	/**
	 * Stores the input path as the log path.
	 * @param logFile	The path will be set to the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
	
	/**
	 * Stores input two strings as names of fastqFiles
	 * @param fastqFile1	The string will be set to the field m_fastqFile1
	 * @param fastqFile2	The string will be set to the field m_fastqFile1
	 */
	public void setFastqFiles(String fastqFile1, String fastqFile2)
	{
		m_fastqFile1 = fastqFile1;
		m_fastqFile2 = fastqFile2;
	}
}
