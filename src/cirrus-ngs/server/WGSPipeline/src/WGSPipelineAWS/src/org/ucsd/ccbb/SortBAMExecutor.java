package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * SortBAMExecutor sorts BAM files.
 * 
 * @author Guorong Xu
 */
public class SortBAMExecutor extends CommandExecutor {
	
	public String m_shellScript = "";
	public String m_inputFile = "";
	public String m_logFile = "";
	
	/**
	 * Assembles a shell command with shell script path, input 
	 * file name and log file path.
	 * This function overrides the function assembleCommands in the file 
	 * CommandExecutor. 
	 * @param filePath	unused in this function
	 * @param fileName	unused in this function
	 * @return a string array includes the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(m_inputFile);
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
	 * Stores the file name as the input file name
	 * @param inputFile		The string will be set to the field m_inputFile
	 */
	public void setInputFile(String inputFile)
	{
		m_inputFile = inputFile;
	}
	
	/**
	 * Stores the input path as the log path.
	 * @param logFile	The path will be set to the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
}
