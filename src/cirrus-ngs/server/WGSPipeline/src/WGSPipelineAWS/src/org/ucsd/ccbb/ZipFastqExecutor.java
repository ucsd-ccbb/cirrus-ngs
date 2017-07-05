package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * ZipFastqExecutor compresses fastq format files.
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class ZipFastqExecutor extends CommandExecutor {
	
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	
	/**
	 * Assembles a shell command with specified shell script path, file name, 
	 * file path and log path.
	 * The function overrides the assembleCommands from the class 
	 * CommandExecutor
	 * @param filePath	a file path that will be added to the command
	 * @param fileName	a file name that will be added to the command
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(filePath + fileName);
		commandArray.add(m_logFile);
		
		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * Stores the input path as the shell script path.
	 * @param fileName	The path will be stored in the field m-shellScript
	 */
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	/**
	 * Stores the input path as the log path.
	 * @param logFile	The path will be stored in the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
}
