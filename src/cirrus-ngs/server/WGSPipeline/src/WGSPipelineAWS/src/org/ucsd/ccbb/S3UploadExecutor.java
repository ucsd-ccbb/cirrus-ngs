package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * S3UploadExecutor uploads files to S3.
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class S3UploadExecutor extends CommandExecutor {
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	
	/**
	 * Assembles a command with specified shell script path, file name, S3 
	 * upload URL and log path.
	 * The function overrides the assembleCommands from the class 
	 * CommandExecutor
	 * @param filePath	unused in this function
	 * @param fileName	a file name that will be added to the command
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(m_fileName);
		commandArray.add(m_S3UploadURL);
		commandArray.add(m_logFile);
		
		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * Stores the input path as the shell script path.
	 * @param fileName	The path will be stored to the field m-shellScript
	 */
	public void setShellFile(String shellScript)
	{
		m_shellScript = shellScript;
	}
	
	/**
	 * Stores the input file name
	 * @param fileName	The file name will be stored in the field m_fileName
	 */
	public void setFileName(String fileName)
	{
		m_fileName = fileName;
	}
	
	/**
	 * Stores the S3 upload URL
	 * @param S3UploadURL	The string will be stored in the field 
	 * 						m_S3UploadURL
	 */
	public void setS3UploadURL(String S3UploadURL)
	{
		m_S3UploadURL = S3UploadURL;
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
