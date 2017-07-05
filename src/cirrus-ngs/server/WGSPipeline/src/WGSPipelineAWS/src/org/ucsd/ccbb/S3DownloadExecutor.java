package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * The class downloads file from S3.
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class S3DownloadExecutor extends CommandExecutor {
	
	/**a string stores the path of the shell script*/
	public String m_shellScript = "";
	/**a string stores the name of a file*/
	public String m_fileName = "";
	/**a string stores the S3 download URL*/
	public String m_s3URL = "";
	/**a string stores the name*/
	public String m_logFile = "";
	
	/**
	 * Assembles a shell command with shell script path, file name, S3 download 
	 * URL and log file path.
	 * @param filePath	the input parameter is not used in this function 
	 * @param fileName	a file name that will be added to the command
	 * @return a string array includes the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(m_fileName);
		commandArray.add(m_s3URL);
		commandArray.add(m_logFile);
		
		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * Stores the input path as the shell script path.
	 * @param fileName	The path will be set to the field m-shellScript
	 */
	public void setShellFile(String shellScript)
	{
		m_shellScript = shellScript;
	}
	
	/**
	 * Stores the input string as the name of a file
	 * @param fileName	The string will be set to the field m_fileName
	 */
	public void setFileName(String fileName)
	{
		m_fileName = fileName;
	}
	
	/**
	 * Stores the input string as the S3 download URL
	 * @param s3URL	The string will be set to the field m_s3URL
	 */
	public void setS3DownloadURL(String s3URL)
	{
		m_s3URL = s3URL;
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
