package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * FastqcExecutor executes fastQC.
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class FastqcExecutor extends CommandExecutor {
	
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the name of the fastq file*/
	public String m_fastqFile = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	
	/**
	 * The function assembles a shell command with shell script path,the fastq
	 * file path and log file path
	 * @param filePath	unused in the function
	 * @param fileName	unused in the function
	 * @return a string array includes the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		
		commandArray.add(m_shellScript);
		commandArray.add(m_fastqFile);
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
	 * Stores the input path as the log path.
	 * @param logFile	The path will be set to the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
	
	/**
	 * Stores the input fastq format file path
	 * @param fastqFile		thr string will be stored in the field m_fastqFile
	 */
	public void setFastqFile(String fastqFile)
	{
		m_fastqFile = fastqFile;
	}
}
