package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * BAM2FastqExecutor assembles a shell command to convert BAM file to fastq file(s).
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class BAM2FastqExecutor extends CommandExecutor {
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the log file path*/
	public String m_logFile = "";
	
	/**
	 * Assembles a shell command with shell script path, file path, 
	 * file name and log file path
	 * The function overrides the function assembleCommands() from class 
	 * CommandExecutor.
	 * @param filePath	a file path that will be added to the command
	 * @param fileName	a file name that will be added to the command
	 * @return a string array of the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		/*create an empty array list of Strings */
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(filePath + fileName);
		commandArray.add(m_logFile);
		
		/*Returns an array containing all the elements of commandArray in proper sequence*/
		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * This function sets the input path as the shell script path.
	 * @param fileName	The path will be set to the field m-shellScript
	 */
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	/**
	 * This function set the input path as the log file path.
	 * @param logFile	The path will be set to the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
}
